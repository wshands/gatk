package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;

import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

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
public final class SVCluster extends MultiVariantWalker {
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String DEFRAG_PADDING_FRACTION_LONG_NAME = "defrag-padding-fraction";
    public static final String CONVERT_INV_LONG_NAME = "convert-inv-to-bnd";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";

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
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    private String variantPrefix = null;

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

    @Argument(
            doc = "Fast mode. Drops hom-ref and no-call genotype fields and emits them as no-calls.",
            fullName = FAST_MODE_LONG_NAME,
            optional = true
    )
    private boolean fastMode = false;

    @Argument(
            doc = "Omit cluster member ID annotations",
            fullName = OMIT_MEMBERS_LONG_NAME,
            optional = true
    )
    private boolean omitMembers = false;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for a breakpoint cluster.",
            optional = true)
    private SVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy = SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName = DEFRAG_PADDING_FRACTION_LONG_NAME,
            doc = "Padding as a fraction of variant length (CNV defragmentation only)",
            optional = true
    )
    private double defragPaddingFraction = CNVDefragmenter.DEFAULT_PADDING_FRACTION;

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    private CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

    @Argument(fullName = JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME,
            doc = "Extend events by this fraction on each side when determining overlap to merge",
            optional = true
    )
    private double defragmentationPadding = CNVDefragmenter.DEFAULT_PADDING_FRACTION;

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
        samples = getSamplesForVariants();

        if (algorithm == CLUSTER_ALGORITHM.DEFRAGMENT_CNV) {
            clusterEngine =  new CNVDefragmenter(dictionary, defragPaddingFraction, clusterParameters.getDepthParameters().getSampleOverlap());
        } else if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVCollapser collapser = new SVCollapser(breakpointSummaryStrategy);
            final LocatableClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ? LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            clusterEngine = new SVClusterEngine<>(dictionary, type, enableCnv, collapser::collapse);
            clusterEngine.setDepthOnlyParams(clusterParameters.getDepthParameters());
            clusterEngine.setMixedParams(clusterParameters.getMixedParameters());
            clusterEngine.setEvidenceParams(clusterParameters.getPESRParameters());
        } else {
            throw new UnsupportedOperationException("Unsupported algorithm: " + algorithm.name());
        }

        writer = createVCFWriter(Paths.get(outputFile));
        writer.writeHeader(createHeader());
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        processClusters();
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
        SVCallRecord call = SVCallRecordUtils.create(variant);
        if (fastMode) {
            final GenotypesContext filteredGenotypes = GenotypesContext.copy(call.getGenotypes().stream()
                    .filter(SVCallRecordUtils::isAltGenotype)
                    .collect(Collectors.toList()));
            call = SVCallRecordUtils.copyCallWithNewGenotypes(call, filteredGenotypes);
        }

        // Flush clusters if we hit the next contig
        if (!call.getContigA().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = call.getContigA();
        }

        // Add to clustering buffer
        if (convertInversions) {
            SVCallRecordUtils.convertInversionsToBreakends(call).forEachOrdered(clusterEngine::add);
        } else {
            clusterEngine.add(call);
        }
    }

    private void processClusters() {
        logger.info("Processing contig " + currentContig + "...");
        List<SVCallRecord> output = clusterEngine.getOutput();
        logger.info("Writing to file...");
        write(output);
        logger.info("Contig " + currentContig + " completed.");
    }

    private void write(final List<SVCallRecord> calls) {
        calls.stream().sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setVCFHeaderVersion(VCFHeaderVersion.VCF4_2);
        header.setSequenceDictionary(dictionary);

        // Copy from inputs
        getHeaderForVariants().getFormatHeaderLines().forEach(header::addMetaDataLine);
        getHeaderForVariants().getInfoHeaderLines().forEach(header::addMetaDataLine);

        // Required info lines
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1, VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine("##INFO=<ID=" + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE + ",Number=.,Type=String,Description=\"Source algorithms\">", header.getVCFHeaderVersion()));
        if (!omitMembers) {
            header.addMetaDataLine(new VCFInfoHeaderLine("##INFO=<ID=" + GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY + ",Number=.,Type=String,Description=\"Cluster variant ids\">", header.getVCFHeaderVersion()));
        }

        // Required format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

        return header;
    }

    public VariantContext buildVariantContext(final SVCallRecord call) {
        final List<Allele> missingAlleles = Arrays.asList(Allele.NO_CALL);
        final GenotypesContext filledGenotypes = SVCallRecordUtils.fillMissingSamplesWithGenotypes(call.getGenotypes(), missingAlleles, samples, Collections.emptyMap());
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsWritten++);
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(), call.getContigB(),
                call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(), call.getAlgorithms(), call.getAlleles(),
                filledGenotypes, call.getAttributes());
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        if (omitMembers) {
            builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        }
        return builder.make();
    }

}
