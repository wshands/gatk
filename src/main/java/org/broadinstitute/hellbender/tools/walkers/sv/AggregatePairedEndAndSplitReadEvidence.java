package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
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
import org.broadinstitute.hellbender.tools.sv.cluster.SVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.SVDeduplicator;
import org.broadinstitute.hellbender.tools.sv.cluster.SVPreprocessingRecordWithEvidenceCollapser;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
    public static final String PE_INNER_WINDOW_LONG_NAME = "pe-inner-window";
    public static final String PE_OUTER_WINDOW_LONG_NAME = "pe-outer-window";
    public static final String SR_WINDOW_LONG_NAME = "sr-window";

    @Argument(
            doc = "Split reads evidence file",
            fullName = SPLIT_READ_LONG_NAME,
            optional = true
    )
    private GATKPath splitReadsFile;

    @Argument(
            doc = "Discordant pairs evidence file",
            fullName = DISCORDANT_PAIRS_LONG_NAME,
            optional = true
    )
    private GATKPath discordantPairsFile;

    @Argument(
            doc = "Tab-delimited table with sample IDs in the first column and expected per-base coverage per sample " +
                    "in the second column. Required if a split read evidence file is provided.",
            fullName = SAMPLE_COVERAGE_LONG_NAME,
            optional = true
    )
    private GATKPath sampleCoverageFile;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Inner discordant pair window size (bp)",
            fullName = PE_INNER_WINDOW_LONG_NAME,
            optional = true
    )
    private int innerWindow = 50;

    @Argument(
            doc = "Outer discordant pair window size (bp)",
            fullName = PE_OUTER_WINDOW_LONG_NAME,
            optional = true
    )
    private int outerWindow = 500;

    @Argument(
            doc = "Split read window size (bp)",
            fullName = SR_WINDOW_LONG_NAME,
            optional = true
    )
    private int splitReadWindow = 200;

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;
    private BreakpointRefiner breakpointRefiner;
    private CachingSVEvidenceAggregator discordantPairCollector;
    private CachingSVEvidenceAggregator startSplitCollector;
    private CachingSVEvidenceAggregator endSplitCollector;
    private SVDeduplicator<SVCallRecordWithEvidence> deduplicator;
    private Map<String,Double> sampleCoverageMap;
    private Set<String> samples;

    private List<SVCallRecordWithEvidence> records;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = new LinkedHashSet<>(getHeaderForVariants().getSampleNamesInOrder());

        if (discordantPairsFile == null && splitReadsFile == null) {
            throw new UserException.BadInput("At least one evidence file must be provided");
        }
        if (discordantPairsFile != null) {
            initializeDiscordantPairDataSource();
            discordantPairCollector = new DiscordantPairEvidenceAggregator(discordantPairSource, dictionary, progressMeter, innerWindow, outerWindow);
        }
        if (splitReadsFile != null) {
            if (sampleCoverageFile == null) {
                throw new UserException.BadInput("Sample coverage table is required when a split read evidence file is present");
            }
            initializeSplitReadEvidenceDataSource();
            loadSampleCoverage();
            startSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, progressMeter, splitReadWindow, true);
            endSplitCollector = new SplitReadEvidenceAggregator(splitReadSource, dictionary, progressMeter, splitReadWindow, false);
            breakpointRefiner = new BreakpointRefiner(sampleCoverageMap, dictionary);
            final SVCollapser<SVCallRecordWithEvidence> collapser = new SVPreprocessingRecordWithEvidenceCollapser(null);
            deduplicator = new SVDeduplicator<>(collapser::collapse, dictionary);
        }
        records = new ArrayList<>();
        writer = createVCFWriter(Paths.get(outputFile));
        writeVCFHeader();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (splitReadSource != null) {
            splitReadSource.close();
        }
        if (discordantPairSource != null) {
            discordantPairSource.close();
        }
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
        records.add(new SVCallRecordWithEvidence(SVCallRecordUtils.create(variant)));
    }

    @Override
    public Object onTraversalSuccess() {
        Stream<SVCallRecordWithEvidence> recordStream = records.stream();

        // Get PE evidence
        if (discordantPairCollector != null) {
            logger.info("Aggregating discordant pair evidence...");
            recordStream = discordantPairCollector.collectEvidence(recordStream.collect(Collectors.toList()));
        }

        // Get SR evidence
        if (startSplitCollector != null && endSplitCollector != null) {
            logger.info("Aggregating start position split-read evidence...");
            recordStream = startSplitCollector.collectEvidence(recordStream.collect(Collectors.toList()));
            logger.info("Aggregating end position split-read evidence...");
            final List<SVCallRecordWithEvidence> endSortedRecords = recordStream
                    .sorted(SVCallRecordUtils.getSVLocatableComparatorByEnds(dictionary))
                    .collect(Collectors.toList());
            recordStream = endSplitCollector.collectEvidence(recordStream.collect(Collectors.toList()));
        }

        // Refine breakpoints
        if (breakpointRefiner != null) {
            recordStream = recordStream.map(r -> breakpointRefiner.refineCall(r));
        }

        // Sort and deduplicate
        recordStream = recordStream.sorted(SVCallRecordUtils.getCallComparator(dictionary));
        if (deduplicator != null) {
            logger.info("Sorting and deduplicating records...");
            recordStream = deduplicator.deduplicateSortedItems(recordStream.collect(Collectors.toList())).stream();
        }

        // Write
        logger.info("Writing to file...");
        write(recordStream);
        return super.onTraversalSuccess();
    }

    private void write(final Stream<SVCallRecordWithEvidence> recordStream) {
        recordStream.map(this::buildVariantContext).forEachOrdered(writer::add);
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
