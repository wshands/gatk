package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Retrieves PE/SR evidence and performs breakpoint refinement
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 *     <li>
 *         Mean depth table
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         SV VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk AggregatePairedEndAndSplitReadEvidence
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Adds PE/SR evidence to structural variant records",
        oneLineSummary = "Adds PE/SR evidence to structural variant records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class AggregatePairedEndAndSplitReadEvidence extends VariantWalker {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";

    @Argument(
            doc = "Split reads evidence file",
            fullName = SPLIT_READ_LONG_NAME
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file",
            fullName = DISCORDANT_PAIRS_LONG_NAME
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first column and expected per-base coverage per sample in the second column",
            fullName = SAMPLE_COVERAGE_LONG_NAME
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private SVDepthOnlyCallDefragmenter defragmenter;
    private BreakpointRefiner breakpointRefiner;
    private PairedEndAndSplitReadEvidenceAggregator evidenceCollector;
    private SVCallRecordDeduplicator<SVCallRecordWithEvidence> deduplicator;
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;
    private int numVariantsWritten = 0;

    private List<SVCallRecord> records;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());
        loadSampleCoverage();
        initializeSplitReadEvidenceDataSource();
        initializeDiscordantPairDataSource();

        defragmenter = new SVDepthOnlyCallDefragmenter(dictionary);
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap, dictionary);
        evidenceCollector = new PairedEndAndSplitReadEvidenceAggregator(splitReadSource, discordantPairSource, dictionary, null);

        final Function<Collection<SVCallRecordWithEvidence>,SVCallRecordWithEvidence> collapser = items -> SVCallRecordUtils.deduplicateWithRawCallAttributeWithEvidence(items, SVCallRecordUtils.ALLELE_COLLAPSER_DIPLOID_NO_CALL);
        deduplicator = new SVCallRecordDeduplicator<>(collapser, dictionary);
        records = new ArrayList<>();

        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    private void initializeSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile.toString(),
                "splitReadsFile",
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void initializeDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile.toString(),
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        final String fileString = sampleCoverageFile.toString();
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(fileString), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(fileString, e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        records.add(SVCallRecordUtils.create(variant));
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Aggregating evidence...");
        final List<SVCallRecordWithEvidence> callsWithEvidence = evidenceCollector.collectEvidence(records);
        logger.info("Filtering and refining breakpoints...");
        final List<SVCallRecordWithEvidence> refinedCalls = callsWithEvidence.stream()
                .map(breakpointRefiner::refineCall)
                .sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .collect(Collectors.toList());
        logger.info("Deduplicating variants...");
        final List<SVCallRecordWithEvidence> finalCalls = deduplicator.deduplicateItems(refinedCalls);
        logger.info("Writing to file...");
        write(finalCalls);

        return super.onTraversalSuccess();
    }

    private void write(final List<SVCallRecordWithEvidence> calls) {
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
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        final Map<String, Object> nonCallAttributes = Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        final List<Allele> missingAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
        final GenotypesContext filledGenotypes = SVCallRecordUtils.fillMissingSamplesWithGenotypes(call.getGenotypes(), missingAlleles, samples, nonCallAttributes);
        final SVCallRecordWithEvidence finalCall = new SVCallRecordWithEvidence(call.getId(), call.getContigA(), call.getPositionA(), call.getStrandA(), call.getContigB(),
                call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(), call.getAlgorithms(),
                filledGenotypes, call.getStartSplitReadSites(), call.getEndSplitReadSites(), call.getDiscordantPairs(),
                call.getCopyNumberDistribution());
        return SVCallRecordUtils.createBuilderWithEvidence(finalCall).make();
    }

}
