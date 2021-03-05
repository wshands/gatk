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
import org.broadinstitute.barclay.argparser.ArgumentCollection;
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

import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME;

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
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String DEFRAG_PADDING_FRACTION_LONG_NAME = "defrag-padding-fraction";
    public static final String CONVERT_INV_LONG_NAME = "convert-inv-to-bnd";
    public static final String ALGORITHM_LONG_NAME = "algorithm";

    enum CLUSTER_ALGORITHM {
        DEFRAGMENT_CNV,
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

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

    @Argument(
            doc = "Enable clustering DEL/DUP variants together as CNVs (does not apply to CNV defragmentation)",
            fullName = ENABLE_CNV_LONG_NAME,
            optional = true
    )
    private boolean enableCnv = false;

    @Argument(
            doc = "Convert inversions to pairs of BNDs",
            fullName = CONVERT_INV_LONG_NAME,
            optional = true
    )
    private boolean convertInversions = false;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for a breakpoint cluster.",
            optional = true)
    private SVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy = SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName = JointGermlineCNVSegmentation.MIN_SAMPLE_NUM_OVERLAP_LONG_NAME,
            doc = "Minimum fraction of common samples for two variants to cluster together (CNV defragmentation only)",
            optional = true
    )
    private double minSampleSetOverlap = CNVDefragmenter.getDefaultSampleOverlap();

    @Argument(fullName = DEFRAG_PADDING_FRACTION_LONG_NAME,
            doc = "Padding as a fraction of variant length (CNV defragmentation only)",
            optional = true
    )
    private double defragPaddingFraction = CNVDefragmenter.getDefaultPaddingFraction();

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    private CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

    @Argument(fullName = JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME,
            doc = "Extend events by this fraction on each side when determining overlap to merge",
            optional = true
    )
    private double defragmentationPadding = CNVDefragmenter.getDefaultPaddingFraction();

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameters = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private SVClusterEngine<SVCallRecord> clusterEngine;
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

        if (algorithm == CLUSTER_ALGORITHM.DEFRAGMENT_CNV) {
            clusterEngine =  new CNVDefragmenter(dictionary, defragPaddingFraction, minSampleSetOverlap);
        } else if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVCollapser<SVCallRecord> collapser = new SVPreprocessingRecordCollapser(breakpointSummaryStrategy);
            final LocatableClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ? LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            clusterEngine = new SVClusterEngine<>(dictionary, type, enableCnv, collapser::collapse);
            clusterEngine.setDepthOnlyParams(clusterParameters.getDepthParameters());
            clusterEngine.setMixedParams(clusterParameters.getMixedParameters());
            clusterEngine.setEvidenceParams(clusterParameters.getPESRParameters());
        } else {
            throw new UnsupportedOperationException("Unsupported algorithm: " + algorithm.name());
        }

        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        if (!clusterEngine.isEmpty()) {
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

        // Add to clustering buffer
        clusterEngine.add(call);
    }

    private void processClusters() {
        logger.info("Processing contig " + currentContig + "...");
        Stream<SVCallRecord> outputStream = clusterEngine.getOutput().stream();
        if (convertInversions) {
            outputStream = outputStream.flatMap(SVCallRecordUtils::convertInversionsToBreakends);
        }
        logger.info("Writing to file...");
        write(outputStream);
        logger.info("Contig " + currentContig + " completed");
    }

    private void write(final Stream<SVCallRecord> calls) {
        calls.sorted(SVCallRecordUtils.getCallComparator(dictionary))
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
