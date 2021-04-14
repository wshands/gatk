package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.lang3.tuple.MutablePair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDB;
import org.broadinstitute.hellbender.utils.samples.Sex;
import org.broadinstitute.hellbender.utils.variant.GATKSVVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Merge GCNV segments VCFs
 *
 * This tool takes in segmented VCFs produced by {@link PostprocessGermlineCNVCalls}.  For single-sample input VCFs,
 * the tool defragments CNV calls by merging adjacent events within the default or specified fractional padding margins.
 * For distinct samples, the tool clusters events that meet the reciprocal overlap requirement. Output is a multi-sample VCF
 * with site allele counts and frequencies (AC, AF).
 *
 * This tool can also be run on multiple multi-sample VCFs, as in a hierarchical gather for large cohorts.  In this case,
 * defragmentation is not preformed because it is assumed that single-sample defragmentation occurred during the initial
 * tool invocation that produced the multi-sample VCF.
 *
 * <h3>Required inputs</h3>
 * <ul>
 *     <li>Reference fasta</li>
 *     <li>Segmented VCFs, single-sample or multi-sample, from {@link PostprocessGermlineCNVCalls}</li>
 *     <li>The interval list used by {@link GermlineCNVCaller} for calling the input VCFs</li>
 *     <li>A pedigree file with an entry giving the sex of each sample</li>
 * </ul>
 *
 * <h3>Output</h3>
 * A single multi-sample VCF with genotypes and copy numbers
 *
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *   gatk JointGermlineCNVSegmentation \
 *         -R reference.fasta
 *         -V sample1.vcf.gz
 *         -V sample2.vcf.gz
 *         --model-call-intervals .filtered.interval_list
 *         --pedigree XandYassignments.ped
 *         -O clustered.vcf.gz
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p><ul>
 * <li>All input samples are required to be present in the pedigree file</li>
 * <li>Quality scores are not output, but can be recalculated based on the updated event boundaries with {@link PostprocessGermlineCNVCalls}</li>
 * <li>Multi-sample input VCFs are assumed to have been generated by this tool and already defragmented</li>
 * <li>This tool only supports mammalian genomes with XX/XY sex determination.</li>
 * </ul></p>
 **/
@BetaFeature
@CommandLineProgramProperties(
        summary = "Gathers single-sample or multi-sample segmented gCNV VCFs, harmonizes breakpoints, and outputs a cohort VCF with genotypes.",
        oneLineSummary = "Combine segmented gCNV VCFs.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class JointGermlineCNVSegmentation extends MultiVariantWalkerGroupedOnStart {

    private SortedSet<String> samples;
    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary dictionary;
    private CNVDefragmenter defragmenter;
    private SVClusterEngine<SVCallRecord> clusterEngine;
    private List<GenomeLoc> callIntervals;
    private String currentContig;
    private SampleDB sampleDB;
    private boolean isMultiSampleInput = false;
    private ReferenceSequenceFile reference;
    private final Set<String> allosomalContigs = new LinkedHashSet<>(Arrays.asList("X","Y","chrX","chrY"));

    class CopyNumberAndEndRecord {
        private MutablePair<Integer, Integer> record;

        public CopyNumberAndEndRecord(final int copyNumber, final int end) {
            record = new MutablePair<>(copyNumber, end);
        }

        public int getCopyNumber() {
            return record.getLeft();
        }

        public int getEndPosition() {
            return record.getRight();
        }
    }

    public static final String MIN_QUALITY_LONG_NAME = "minimum-qs-score";
    public static final String MIN_SAMPLE_NUM_OVERLAP_LONG_NAME = "min-sample-set-fraction-overlap";
    public static final String DEFRAGMENTATION_PADDING_LONG_NAME = "defragmentation-padding-fraction";
    public static final String CLUSTERING_INTERVAL_OVERLAP_LONG_NAME = "clustering-interval-overlap";
    public static final String CLUSTERING_BREAKEND_WINDOW_LONG_NAME = "clustering-breakend-window";
    public static final String MODEL_CALL_INTERVALS_LONG_NAME = "model-call-intervals";
    public static final String BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME = "breakpoint-summary-strategy";

    @Argument(fullName = MIN_QUALITY_LONG_NAME, doc = "Minimum QS score to combine a variant segment", optional = true)
    private int minQS = 20;

    @Argument(fullName = MIN_SAMPLE_NUM_OVERLAP_LONG_NAME, doc = "Minimum fraction of common samples for two variants to cluster together", optional = true)
    private double minSampleSetOverlap = CNVDefragmenter.getDefaultSampleOverlap();

    @Argument(fullName = DEFRAGMENTATION_PADDING_LONG_NAME, doc = "Extend events by this fraction on each side when determining overlap to merge", optional = true)
    private double defragmentationPadding = CNVDefragmenter.getDefaultPaddingFraction();

    @Argument(fullName = CLUSTERING_INTERVAL_OVERLAP_LONG_NAME,
            doc="Minimum interval reciprocal overlap for clustering", optional=true)
    public double clusterIntervalOverlap = SVClusterEngine.DEFAULT_DEPTH_ONLY_PARAMS.getReciprocalOverlap();

    @Argument(fullName = CLUSTERING_BREAKEND_WINDOW_LONG_NAME,
            doc="Cluster events whose endpoints are within this distance of each other", optional=true)
    public int clusterWindow = SVClusterEngine.DEFAULT_DEPTH_ONLY_PARAMS.getWindow();

    @Argument(fullName = MODEL_CALL_INTERVALS_LONG_NAME, doc = "gCNV model intervals created with the FilterIntervals tool.")
    private GATKPath modelCallIntervalList = null;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME, doc = "Strategy to use for choosing a representative value for a breakpoint cluster.", optional = true)
    private SVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy = SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined output VCF")
    private GATKPath outputFile;

    @Argument(
            doc = "Reference copy-number on autosomal intervals.",
            fullName = PostprocessGermlineCNVCalls.AUTOSOMAL_REF_COPY_NUMBER_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int refAutosomalCopyNumber = 2;

    /**
     * See https://software.broadinstitute.org/gatk/documentation/article.php?id=7696 for more details on the PED
     * format. Note that each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to
     * tell the GATK PED parser that the corresponding fields are missing from the ped file.
     *
     */
    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for samples")
    private GATKPath pedigreeFile = null;

    @Override
    public boolean doDictionaryCrossValidation() {
        return false;
    }

    //require a reference to do dictionary validation since there may be too many samples for cross-validating
    @Override
    public boolean requiresReference() {
        return true;
    }

    // Cannot require sample overlap when clustering across samples
    private static final double CLUSTER_SAMPLE_OVERLAP_FRACTION = 0;

    @Override
    public void onTraversalStart() {
        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());

        //strict validation will verify that all samples are in the pedigree
        sampleDB = SampleDB.createSampleDBFromPedigreeAndDataSources(pedigreeFile, getSamplesForVariants(), PedigreeValidationType.STRICT);

        dictionary = getBestAvailableSequenceDictionary();
        //dictionary will not be null because this tool requiresReference()

        final GenomeLocParser parser = new GenomeLocParser(this.dictionary);
        setIntervals(parser);

        if (callIntervals == null) {
            defragmenter = new CNVDefragmenter(dictionary, defragmentationPadding, minSampleSetOverlap);
        } else {
            defragmenter = new BinnedCNVDefragmenter(dictionary, defragmentationPadding, minSampleSetOverlap, callIntervals);
        }
        clusterEngine = new SVClusterEngine<>(dictionary, LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE,
                true, (new SVCollapser(breakpointSummaryStrategy))::collapse);
        final SVClusterEngine.DepthClusteringParameters clusterArgs = new SVClusterEngine.DepthClusteringParameters(clusterIntervalOverlap, clusterWindow, CLUSTER_SAMPLE_OVERLAP_FRACTION);
        clusterEngine.setDepthOnlyParams(clusterArgs);

        vcfWriter = getVCFWriter();

        if (getSamplesForVariants().size() != 1) {
            logger.warn("Multi-sample VCFs found, which are assumed to be pre-clustered. Skipping defragmentation.");
            isMultiSampleInput = true;
        } else {
            isMultiSampleInput = false;
        }
    }

    /**
     * If model intervals are supplied, subset to the requested traversal intervals
     * @param parser    needed to merge intervals if necessary
     */
    private void setIntervals(final GenomeLocParser parser) {
        if (modelCallIntervalList != null) {
           final List<GenomeLoc> inputCoverageIntervals = IntervalUtils.featureFileToIntervals(parser, modelCallIntervalList.getURIString());
            final List<GenomeLoc> inputTraversalIntervals = IntervalUtils.genomeLocsFromLocatables(parser,getTraversalIntervals());
            callIntervals = IntervalUtils.mergeListsBySetOperator(inputCoverageIntervals, inputTraversalIntervals, IntervalSetRule.INTERSECTION);
        }
    }

    private VariantContextWriter getVCFWriter() {
        samples = getSamplesForVariants();

        final VCFHeader inputVCFHeader = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        headerLines.add(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));

        VariantContextWriter writer = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = new IndexedSampleList(samples).asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        writer.writeHeader(vcfHeader);

        return writer;
    }

    /**
     * @param variantContexts  VariantContexts from driving variants with matching start position
     *                         NOTE: This will never be empty
     * @param referenceContext ReferenceContext object covering the reference of the longest spanning VariantContext
     * @param readsContexts
     */
    @Override
    public void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        if (currentContig == null) {
            currentContig = variantContexts.get(0).getContig(); //variantContexts should have identical start, so choose 0th arbitrarily
        } else if (!variantContexts.get(0).getContig().equals(currentContig)) {
            processClusters();
            currentContig = variantContexts.get(0).getContig();
        }
        for (final VariantContext vc : variantContexts) {
            final SVCallRecord record = SVCallRecordUtils.createDepthOnlyFromGCNVWithOriginalGenotypes(vc, minQS);
            if (record != null) {
                if (!isMultiSampleInput) {
                    defragmenter.add(record);
                } else {
                    clusterEngine.add(record);
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        processClusters();
        return null;
    }

    private void processClusters() {
        if (!defragmenter.isEmpty()) {
            final List<SVCallRecord> defragmentedCalls = defragmenter.getOutput();
            defragmentedCalls.stream().forEachOrdered(clusterEngine::add);
        }
        //Jack and Isaac cluster first and then defragment
        final List<SVCallRecord> clusteredCalls = clusterEngine.getOutput();
        write(clusteredCalls);
    }

    private void write(final List<SVCallRecord> calls) {
        final List<VariantContext> sortedCalls = calls.stream()
                .sorted(Comparator.comparing(c -> new SimpleInterval(c.getContigA(), c.getPositionA(), c.getPositionB()), //VCs have to be sorted by end as well
                        IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(record -> buildVariantContext(record, reference))
                .collect(Collectors.toList());
        final Iterator<VariantContext> it = sortedCalls.iterator();
        ArrayList<VariantContext> overlappingVCs = new ArrayList<>(calls.size());
        if (!it.hasNext()) {
            return;
        }
        int clusterEnd = -1;
        String clusterContig = null;
        //gather groups of overlapping VCs and update the genotype copy numbers appropriately
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            if ((clusterEnd == -1 || curr.getStart() < clusterEnd) && (clusterContig == null || curr.getContig().equals(clusterContig))) {
                overlappingVCs.add(curr);
                if (curr.getEnd() > clusterEnd) {
                    clusterEnd = curr.getEnd();
                }
                if (clusterContig == null) {
                    clusterContig = curr.getContig();
                }
            } else {
                final List<VariantContext> resolvedVCs = resolveVariantContexts(allosomalContigs, refAutosomalCopyNumber, sampleDB, samples, overlappingVCs);
                resolvedVCs.forEach(vcfWriter::add);
                overlappingVCs = new ArrayList<>();
                overlappingVCs.add(curr);
                clusterEnd = curr.getEnd();
                clusterContig = curr.getContig();
            }
        }
        //write out the last set of overlapping VCs
        final List<VariantContext> resolvedVCs = resolveVariantContexts(allosomalContigs, refAutosomalCopyNumber, sampleDB, samples, overlappingVCs);
        resolvedVCs.forEach(vcfWriter::add);
    }

    /**
     * Correct genotype calls for overlapping variant contexts
     * Note that we assume that a sample will not occur twice with the same copy number because it should have been defragmented
     * @param allosomalContigs    names of allosomal contigs (e.g. X, Y, chrX, chrY)
     * @param refAutosomalCopyNumber   reference copy number for autosomes
     * @param sampleDB data structure containing sample sex assignments
     * @param samples   full set of samples to be output
     * @param overlappingVCs    overlapping variant contexts with correct genotypes, but maybe not copy numbers
     * @return a list of VariantContexts with genotypes corrected for called alleles and copy number
     */
    @VisibleForTesting
    protected List<VariantContext> resolveVariantContexts(final Set<String> allosomalContigs, final int refAutosomalCopyNumber,
                                                          final SampleDB sampleDB, final SortedSet<String> samples,
                                                          final List<VariantContext> overlappingVCs) {
        Utils.nonNull(overlappingVCs);
        final List<VariantContext> resolvedVCs = new ArrayList<>(overlappingVCs.size());
        final Iterator<VariantContext> it = overlappingVCs.iterator();

        //sampleName, copyNumber, endPos -- it's safe to just use position because if the VCs overlap then they must be on the same contig
        final Map<String, CopyNumberAndEndRecord> sampleCopyNumbers = new LinkedHashMap<>(SVUtils.hashMapCapacity(overlappingVCs.size()));
        while (it.hasNext()) {
            final VariantContext curr = it.next();
            resolvedVCs.add(updateGenotypes(allosomalContigs, refAutosomalCopyNumber, sampleDB, samples, curr, sampleCopyNumbers));
            //update copy number table for subsequent VCs using variant genotypes from input VCs
            for (final Genotype g : curr.getGenotypes()) {
                if (g.hasAnyAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                    sampleCopyNumbers.put(g.getSampleName(),
                            new CopyNumberAndEndRecord(
                                    VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, refAutosomalCopyNumber),
                                    curr.getAttributeAsInt(VCFConstants.END_KEY, curr.getStart())));
                }
            }
        }
        return resolvedVCs;
    }

    /**
     * For CNVs, i.e. DEL and DUP alts only
     * @param allosomalContigs  names of allosomal contigs (e.g. X, Y, chrX, chrY)
     * @param refAutosomalCopyNumber    reference copy number for autosomes
     * @param sampleDB  data structure containing sample sex assignments
     * @param samples   full set of samples to be output
     * @param vc VariantContext with just variant samples
     * @param sampleCopyNumbers may be modified to remove terminated variants
     * @return new VariantContext with AC and AF
     */
    @VisibleForTesting
    protected static VariantContext updateGenotypes(final Set<String> allosomalContigs, final int refAutosomalCopyNumber,
                                                    final SampleDB sampleDB, final SortedSet<String> samples,
                                                    final VariantContext vc, final Map<String, CopyNumberAndEndRecord> sampleCopyNumbers) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final List<Genotype> newGenotypes = new ArrayList<>();
        final Allele vcRefAllele = vc.getReference();
        final Map<Allele, Long> alleleCountMap = new LinkedHashMap<>(2);
        if (vc.getAlternateAlleles().stream().filter(a -> !a.equals(GATKSVVCFConstants.DEL_ALLELE)).filter(a -> !a.equals(GATKSVVCFConstants.DUP_ALLELE)).count() > 0) {
            throw new IllegalArgumentException("At site " + vc.getContig() + ":" + vc.getStart() + " variant context contains alternate alleles in addition to CNV <DEL> and <DUP> alleles: " + vc.getAlternateAlleles());
        }
        alleleCountMap.put(GATKSVVCFConstants.DEL_ALLELE, 0L);
        alleleCountMap.put(GATKSVVCFConstants.DUP_ALLELE, 0L);
        int alleleNumber = 0;
        for (final String sample : samples) {
            final Genotype g = vc.getGenotype(sample); //may be null
            final GenotypeBuilder genotypeBuilder = g == null? new GenotypeBuilder(sample) : new GenotypeBuilder(g);  //use g to set alleles
            final List<Allele> alleles;

            //get proper ploidy for autosomes and allosomes (sex chromosomes)
            final int samplePloidy = getSamplePloidy(allosomalContigs, refAutosomalCopyNumber, sampleDB, sample, vc.getContig(), g);
            alleleNumber += samplePloidy;

            //"square off" the genotype matrix by filling in missing (non-variant) samples with homRef calls with reference copy number
            if (!sampleCopyNumbers.containsKey(sample) && !vc.hasGenotype(sample)) {
                genotypeBuilder.alleles(GATKVariantContextUtils.makePloidyLengthAlleleList(samplePloidy, vcRefAllele));
                genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, samplePloidy);
                newGenotypes.add(genotypeBuilder.make());
            //handle variant genotypes and reference genotypes with non-reference copy number due to overlapping events
            } else {
                //determine sample copy number from VC genotype or sampleCopyNumbers map or default to reference, i.e. samplePloidy
                final int copyNumber;
                if (sampleCopyNumbers.containsKey(sample) && sampleCopyNumbers.get(sample).getEndPosition() > vc.getStart()) {
                    copyNumber = sampleCopyNumbers.get(sample).getCopyNumber();
                    alleles = GATKVariantContextUtils.makePloidyLengthAlleleList(samplePloidy, vcRefAllele);
                } else if (g != null) {
                    copyNumber = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, samplePloidy);
                    if (samplePloidy == g.getPloidy()) {
                        alleles = g.getAlleles();
                    } else {
                        alleles = GATKSVVariantContextUtils.makeGenotypeAllelesFromCopyNumber(copyNumber, samplePloidy, vcRefAllele);
                    }
                } else {
                    copyNumber = samplePloidy;
                    alleles = GATKSVVariantContextUtils.makeGenotypeAllelesFromCopyNumber(copyNumber, samplePloidy, vcRefAllele);
                }
                genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);

                genotypeBuilder.alleles(alleles);
                newGenotypes.add(genotypeBuilder.make());

                //use called alleles to update AC
                //check for genotype in VC because we don't want to count overlapping events (in sampleCopyNumbers map) towards AC
                if (vc.hasGenotype(sample)) {
                    if (alleles.contains(GATKSVVCFConstants.DEL_ALLELE)) {
                        final Long count = alleleCountMap.get(GATKSVVCFConstants.DEL_ALLELE);
                        alleleCountMap.put(GATKSVVCFConstants.DEL_ALLELE, count + alleles.stream().filter(Allele::isNonReference).count());
                    } else if (copyNumber > samplePloidy) {
                        final Long count = alleleCountMap.get(GATKSVVCFConstants.DUP_ALLELE);
                        alleleCountMap.put(GATKSVVCFConstants.DUP_ALLELE, count + 1); //best we can do for dupes is carrier frequency
                    }
                }
            }
        }
        builder.genotypes(newGenotypes);

        if (alleleNumber > 0) {
            if (vc.getAlternateAlleles().size() == 1) {
                final long AC = alleleCountMap.get(vc.getAlternateAllele(0));
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, AC)
                        .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, (double)AC / alleleNumber)
                        .attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
            } else {
                final List<Long> alleleCounts = new ArrayList<>(vc.getNAlleles());
                final List<Double> alleleFreqs = new ArrayList<>(vc.getNAlleles());
                for (final Allele a : builder.getAlleles()) {
                    if (a.isReference()) {
                        continue;
                    }
                    alleleCounts.add(alleleCountMap.get(a));
                    alleleFreqs.add(Double.valueOf(alleleCountMap.get(a)));
                }
                builder.attribute(VCFConstants.ALLELE_COUNT_KEY, alleleCounts)
                    .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs)
                    .attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
            }
        }
        return builder.make();
    }

    /**
     * Get the sample ploidy for a given contig based on the input pedigree (if available)
     * defaults to {@param refAutosomalCopyNumber} for autosomes
     * @param allosomalContigs  names of allosomal contigs (e.g. X, Y, chrX, chrY)
     * @param refAutosomalCopyNumber    reference copy number for autosomes
     * @param sampleDB  data structure containing sample sex assignments
     * @param sampleName    current sample of interest
     * @param contig    current contig of interest
     * @param g may be null
     * @return
     */
    @VisibleForTesting
    protected static int getSamplePloidy(final Set<String> allosomalContigs, final int refAutosomalCopyNumber,
                                       final SampleDB sampleDB, final String sampleName, final String contig, final Genotype g) {
        if (!allosomalContigs.contains(contig)) {
            return refAutosomalCopyNumber;
        }
        if (sampleDB == null || sampleDB.getSample(sampleName) == null) {
            if (g != null) {
                return g.getPloidy();
            } else {
                throw new IllegalStateException("Sample " + sampleName + " is missing from the pedigree");
            }
        } else {
            final Sex sampleSex = sampleDB.getSample(sampleName).getSex();
            if (contig.equals("X") || contig.equals("chrX")) {
                if (sampleSex.equals(Sex.FEMALE)) {
                    return 2;
                } else if (sampleSex.equals(Sex.MALE)) {
                    return 1;
                } else { //UNKNOWN
                    return 1;
                }
            } else if (contig.equals("Y") || contig.equals("chrY")) {
                if (sampleSex.equals(Sex.FEMALE)) {
                    return 0;
                } else if (sampleSex.equals(Sex.MALE)) {
                    return 1;
                } else { //UNKNOWN
                    return 1;
                }
            } else {
                throw new IllegalArgumentException("Encountered unknown allosomal contig: " + contig + ". This tool only " +
                        "supports mammalian genomes with XX/XY sex determination.");
            }
        }
    }

    @VisibleForTesting
    protected static VariantContext buildVariantContext(final SVCallRecord call, final ReferenceSequenceFile reference) {
        Utils.nonNull(call);
        Utils.nonNull(reference);
        final boolean isCNV = call.getType().equals(StructuralVariantType.CNV);
        final List<Allele> outputAlleles = new ArrayList<>(3);  //max is ref, del, dup
        final Allele refAllele = Allele.create(ReferenceUtils.getRefBaseAtPosition(reference, call.getContigA(), call.getPositionA()), true);
        outputAlleles.add(refAllele);
        if (!isCNV) {
            outputAlleles.add(Allele.create("<" + call.getType().name() + ">", false));
        } else {
            outputAlleles.add(GATKSVVCFConstants.DEL_ALLELE);
            outputAlleles.add(GATKSVVCFConstants.DUP_ALLELE);
        }

        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContigA(), call.getPositionA(), call.getPositionB(),
                outputAlleles);
        builder.attribute(VCFConstants.END_KEY, call.getPositionB());
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        if (isCNV) {
            builder.attribute(VCFConstants.SVTYPE, "MCNV");  //MCNV for compatibility with svtk annotate
        } else {
            builder.attribute(VCFConstants.SVTYPE, call.getType());
        }
        final List<Genotype> genotypes = new ArrayList<>(call.getGenotypes().size());
        for (final Genotype g : call.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(g);
            //update reference alleles
            final List<Allele> newGenotypeAlleles = new ArrayList<>(g.getAlleles().size());
            for (final Allele a : g.getAlleles()) {
                if (a.isReference()) {
                    newGenotypeAlleles.add(refAllele);
                } else {
                    newGenotypeAlleles.add(a);
                }
            }
            genotypeBuilder.alleles(newGenotypeAlleles);
             if (g.hasAnyAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
            }
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder.make();
    }

    @Override
    public void closeTool(){
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
