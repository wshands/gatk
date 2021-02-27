package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;

import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

/**
 * Clusters SVs with similar breakpoints based on coordinates and supporting evidence.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Unclustered structural variants from
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Clusters structural variants",
        oneLineSummary = "Clusters structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVCluster extends VariantWalker {
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Prefix for variant IDs",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    private String variantPrefix = "SV_x";

    @Argument(fullName = JointGermlineCNVSegmentation.MIN_SAMPLE_NUM_OVERLAP_LONG_NAME,
            doc = "Minimum fraction of common samples for two variants to cluster together",
            optional = true
    )
    private double minSampleSetOverlap = CNVDefragmenter.getDefaultSampleOverlap();

    @Argument(fullName = JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME,
            doc = "Extend events by this fraction on each side when determining overlap to merge",
            optional = true
    )
    private double defragmentationPadding = CNVDefragmenter.getDefaultPaddingFraction();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private CNVDefragmenter defragmenter;
    private List<SVCallRecord> nonDepthRawCallsBuffer;
    private SVClusterEngine<SVCallRecord> singleLinkageEngine;
    private SVClusterEngine<SVCallRecord> maxCliqueEngine;
    private Set<String> samples;
    private String currentContig;
    private int numVariantsWritten = 0;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());

        defragmenter = new CNVDefragmenter(dictionary, 0.5, 0.9);
        nonDepthRawCallsBuffer = new ArrayList<>();

        final SVCollapser<SVCallRecord> collapser = new SVPreprocessingRecordCollapser(SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END);
        singleLinkageEngine = new SVClusterEngine<>(dictionary, LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, false, collapser::collapse);
        singleLinkageEngine.setDepthOnlyParams(new SVClusterEngine.DepthClusteringParameters(0.9, 0, 0));
        singleLinkageEngine.setMixedParams(new SVClusterEngine.MixedClusteringParameters(0.9, 50, 50));
        singleLinkageEngine.setEvidenceParams(new SVClusterEngine.EvidenceClusteringParameters(0.9, 50, 50));

        maxCliqueEngine = new SVClusterEngine<>(dictionary, LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE, false, collapser::collapse);

        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        if (!defragmenter.isEmpty() || !nonDepthRawCallsBuffer.isEmpty()) {
            processClusters();
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord originalCall = SVCallRecordUtils.create(variant);
        final ArrayList<Genotype> filteredGenotypeList = new ArrayList<>(originalCall.getGenotypes().size());
        originalCall.getGenotypes().stream()
                .filter(SVCallRecord::isRawCall)
                .forEach(filteredGenotypeList::add);
        final GenotypesContext filteredGenotypes = GenotypesContext.create(filteredGenotypeList);
        final SVCallRecord call = SVCallRecordUtils.copyCallWithNewGenotypes(originalCall, filteredGenotypes);

        // Flush clusters if we hit the next contig
        if (!call.getContigA().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = call.getContigA();
        }

        // Add to clustering buffers
        if (CNVDefragmenter.isDepthOnlyCall(call)) {
            defragmenter.add(call);
        } else {
            nonDepthRawCallsBuffer.add(call);
        }
    }

    private void processClusters() {
        logger.info("Processing contig " + currentContig + "...");
        final Stream<SVCallRecord> defragmentedStream = defragmenter.getOutput().stream();
        final Stream<SVCallRecord> nonDepthStream = nonDepthRawCallsBuffer.stream()
                .flatMap(SVCallRecordUtils::convertInversionsToBreakends);
        //Combine and sort depth and non-depth calls because they must be added in dictionary order
        Stream.concat(defragmentedStream, nonDepthStream)
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .forEachOrdered(singleLinkageEngine::add);
        nonDepthRawCallsBuffer.clear();
        logger.info("Clustering with single-linkage...");
        final List<SVCallRecord> singleLinkageCalls = singleLinkageEngine.getOutput();

        logger.info("Clustering with maximal clique...");
        //singleLinkageCalls.stream().forEachOrdered(maxCliqueEngine::add); TODO
        //final List<SVCallRecord> clusteredCalls = maxCliqueEngine.getOutput();

        logger.info("Writing to file...");
        write(singleLinkageCalls);
        logger.info("Contig " + currentContig + " completed");
    }

    private void write(final List<SVCallRecord> calls) {
        calls.stream()
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setSequenceDictionary(dictionary);
        for (final VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            header.addMetaDataLine(line);
        }
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecord call) {
        final Map<String, Object> nonCallAttributes = Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        final List<Allele> missingAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
        final GenotypesContext filledGenotypes = SVCallRecordUtils.fillMissingSamplesWithGenotypes(call.getGenotypes(), missingAlleles, samples, nonCallAttributes);
        final String newId = String.format("%s%08x", variantPrefix, numVariantsWritten++);
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(), call.getContigB(),
                call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(), call.getAlgorithms(),
                filledGenotypes);
        return SVCallRecordUtils.getVariantBuilder(finalCall).make();
    }

}
