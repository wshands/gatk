package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.Permutation;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Condense homRef blocks in a single-sample GVCF
 *
 * <p>
 * ReblockGVCF compresses a GVCF by merging hom-ref blocks that were produced using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller according to new GQ band parameters.  A joint callset produced with GVCFs reprocessed by ReblockGVCF will have
 * lower precision for hom-ref genotype qualities at variant sites, but the input data footprint can be greatly reduced
 * if the default GQ band parameters are used.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A HaplotypeCaller-produced GVCF to reblock
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A smaller GVCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk ReblockGVCF \
 *   -R reference.fasta \
 *   -V sample1.g.vcf \
 *   -O sample1.reblocked.g.vcf
 * </pre>
 *
 * Invocation as for use with GnarlyGenotyper in the "Biggest Practices"
 * <pre>
 *  gatk ReblockGVCF \
 *    -R reference.fasta \
 *    -V sample1.g.vcf \
 *    -drop-low-quals \
 *    -rgq-threshold 10 \
 *    -do-qual-approx \
 *    -O sample1.reblocked.g.vcf
 *  * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only single-sample GVCF files produced by HaplotypeCaller can be used as input for this tool.</p>
 *
 */
@BetaFeature
@CommandLineProgramProperties(summary = "Compress a single-sample GVCF from HaplotypeCaller by merging homRef blocks using new GQ band parameters",
        oneLineSummary = "Condenses homRef blocks in a single-sample GVCF",
        programGroup = OtherProgramGroup.class,
        omitFromCommandLine = true)
@DocumentedFeature
public final class ReblockGVCF extends MultiVariantWalker {

    private static final OneShotLogger logger = new OneShotLogger(ReblockGVCF.class);

    private final static int PLOIDY_TWO = 2;  //assume diploid genotypes

    private int bufferEnd = 0;
    private int vcfOutputEnd = 0;

    public static final String DROP_LOW_QUALS_ARG_NAME = "drop-low-quals";
    public static final String RGQ_THRESHOLD_LONG_NAME = "rgq-threshold-to-no-call";
    public static final String RGQ_THRESHOLD_SHORT_NAME = "rgq-threshold";
    public static final String KEEP_ALL_ALTS_ARG_NAME = "keep-all-alts";
    public static final String QUAL_APPROX_LONG_NAME = "do-qual-score-approximation";
    public static final String QUAL_APPROX_SHORT_NAME = "do-qual-approx";
    public static final String ALLOW_MISSING_LONG_NAME = "allow-missing-hom-ref-data";
    public static final String POSTERIORS_KEY_LONG_NAME = "genotype-posteriors-key";

    private static final Comparator<? super VariantContextBuilder> VCB_COMPARATOR = Comparator.comparingLong(VariantContextBuilder::getStart);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    private GATKPath outputFile;

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Output the band lower bound for each GQ block instead of the min GQ -- for better compression
     */
    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.OUTPUT_BLOCK_LOWER_BOUNDS, doc = "Output the band lower bound for each GQ block regardless of the data it represents", optional = true)
    private boolean floorBlocks = false;

    @Advanced
    @Argument(fullName=HaplotypeCallerArgumentCollection.GQ_BAND_LONG_NAME, shortName=HaplotypeCallerArgumentCollection.GQ_BAND_SHORT_NAME,
            doc="Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)")
    public List<Integer> GVCFGQBands = new ArrayList<>();
    {
        GVCFGQBands.add(20); GVCFGQBands.add(100);
    }

    @Advanced
    @Argument(fullName=DROP_LOW_QUALS_ARG_NAME, shortName=DROP_LOW_QUALS_ARG_NAME, doc="Exclude variants and homRef blocks that are GQ0 from the reblocked GVCF to save space; drop low quality/uncalled alleles")
    protected boolean dropLowQuals = false;

    @Advanced
    @Argument(fullName=RGQ_THRESHOLD_LONG_NAME, shortName=RGQ_THRESHOLD_SHORT_NAME, doc="Reference genotype quality (PL[0]) value below which variant sites will be converted to GQ0 homRef calls")
    protected double rgqThreshold = 0.0;

    @Advanced
    @Argument(fullName=QUAL_APPROX_LONG_NAME, shortName=QUAL_APPROX_SHORT_NAME, doc="Add necessary INFO field annotation to perform QUAL approximation downstream; required for GnarlyGenotyper")
    protected boolean doQualApprox = false;

    @Advanced
    @Argument(fullName=ALLOW_MISSING_LONG_NAME, doc="Fill in homozygous reference genotypes with no PLs and no GQ with PL=[0,0,0].  Necessary for input from Regeneron's WeCall variant caller.")
    protected boolean allowMissingHomRefData = false;

    @Advanced
    @Argument(fullName=KEEP_ALL_ALTS_ARG_NAME, doc="Keep all ALT alleles and full PL array for most accurate GQs")
    protected boolean keepAllAlts = false;

    @Advanced
    @Argument(fullName=POSTERIORS_KEY_LONG_NAME, doc="INFO field key corresponding to the posterior genotype probabilities", optional = true)
    protected String posteriorsKey = null;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;
    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = Arrays.asList(GVCFWriter.GVCF_BLOCK, GATKVCFConstants.HAPLOTYPE_SCORE_KEY,
            GATKVCFConstants.INBREEDING_COEFFICIENT_KEY, GATKVCFConstants.MLE_ALLELE_COUNT_KEY,
            GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, GATKVCFConstants.EXCESS_HET_KEY, GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY,
            GATKVCFConstants.DOWNSAMPLED_KEY);

    private final List<VariantContextBuilder> homRefBlockBuffer = new ArrayList<>(10);  //10 is a generous estimate for the number of overlapping deletions
    private String currentContig;
    private VariantContextWriter vcfWriter;
    private CachingIndexedFastaSequenceFile referenceReader;

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class);
    }

    @Override
    public boolean requiresReference() {return true;}

    @Override
    public void onTraversalStart() {
        if (getSamplesForVariants().size() != 1) {
            throw new UserException.BadInput("ReblockGVCF can take multiple input GVCFs, but they must be "
                    + "non-overlapping shards from the same sample.  Found samples " + getSamplesForVariants());
        }

        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeaders);
        // Remove GCVFBlocks, legacy headers, and annotations that aren't informative for single samples
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCFWriter.GVCF_BLOCK) ||
                (vcfHeaderLine.getKey().equals("INFO")) && ((VCFInfoHeaderLine)vcfHeaderLine).getID().equals(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ||  //remove old (maybe wrong type) and add new with deprecated note
                (vcfHeaderLine.getKey().equals("INFO")) && infoFieldAnnotationKeyNamesToRemove.contains(((VCFInfoHeaderLine)vcfHeaderLine).getID()));

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        genotypingEngine = createGenotypingEngine(new IndexedSampleList(inputHeader.getGenotypeSamples()));
        createAnnotationEngine();

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED));  //NOTE: this is deprecated, but keep until we reprocess all GVCFs
        if (inputHeader.hasInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED));
        }

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        final VariantContextWriter writer = createVCFWriter(outputFile);

        try {
            vcfWriter = new GVCFWriter(writer, new ArrayList<>(GVCFGQBands), floorBlocks);
        } catch ( final IllegalArgumentException e ) {
            throw new IllegalArgumentException("GQBands are malformed: " + e.getMessage(), e);
        }
        vcfWriter.writeHeader(new VCFHeader(headerLines, getSamplesForVariants()));  //don't get samples from header -- multi-variant inputHeader doens't have sample names

        if (genotypeArgs.samplePloidy != PLOIDY_TWO) {
            throw new UserException.BadInput("The -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this tool maintains input ploidy for each genotype");
        }

         referenceReader = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
    }

    private HaplotypeCallerGenotypingEngine createGenotypingEngine(final SampleList samples) {
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        // create the genotyping engine
        hcArgs.standardArgs.outputMode = OutputMode.EMIT_ALL_CONFIDENT_SITES;  //use confident vs. active mode so we can drop low quality variants
        hcArgs.standardArgs.annotateAllSitesWithPLs = true;
        hcArgs.standardArgs.genotypeArgs = genotypeArgs.clone();
        hcArgs.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;   //this is important to force emission of all alleles at a multiallelic site
        hcArgs.standardArgs.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING = dropLowQuals ? genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING : 0.0;
        return new HaplotypeCallerGenotypingEngine(hcArgs, samples, true, false);

    }

    @VisibleForTesting
    protected void createAnnotationEngine() {
        annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), dbsnp.dbsnp, Collections.emptyList(), false, false);
    }

    // get VariantContexts from input gVCFs and regenotype
    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        if (currentContig == null) {
            currentContig = variant.getContig(); //variantContexts should have identical start, so choose 0th arbitrarily
        } else if (!variant.getContig().equals(currentContig)) {
            flushRefBlockBuffer();
            currentContig = variant.getContig();
            vcfOutputEnd = 0;
        }
        final VariantContext newVC;
        try {
            newVC = regenotypeVC(variant);
        } catch (final Exception e) {
            throw new GATKException("Exception thrown at " + variant.getContig() + ":" + variant.getStart() + " " + variant.toString(), e);
        }
        if (newVC != null) {
            try {
                vcfWriter.add(newVC);
                if (newVC.getEnd() > vcfOutputEnd) {
                    vcfOutputEnd = newVC.getEnd();
                }
            } catch (Exception e) {
                throw new GATKException("Exception thrown at " + newVC.getContig() + ":" + newVC.getStart() + " " + newVC.toString(), e);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        flushRefBlockBuffer();
        return null;
    }

    /**
     * Re-genotype (and re-annotate) a VariantContext
     * Note that the GVCF write takes care of the actual homRef block merging based on {@code GVCFGQBands}
     *
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic (may have been added to ref block buffer)
     */
    private VariantContext regenotypeVC(final VariantContext originalVC) {
        VariantContext result = originalVC;

        //Pass back ref-conf homRef sites/blocks to be combined by the GVCFWriter
        if (isHomRefBlock(result)) {
            if (result.getEnd() <= vcfOutputEnd) {
                return null;
            }
            final VariantContext filtered = processHomRefBlock(result);
            if (filtered != null) {
                updateHomRefBlockBuffer(filtered);
            }
            return null;
        }

        //Use the genotyping engine to do the QUAL thresholding if we're dropping low qual sites
        //don't need to calculate quals for sites with no data whatsoever or sites already genotyped homRef,
        //but if STAND_CALL_CONF > 0 we need to drop low quality alleles and regenotype
        //Note that spanning deletion star alleles will be considered low quality
        if (dropLowQuals && result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 0 && !isMonomorphicCallWithAlts(result)) {
            final VariantContext regenotyped = genotypingEngine.calculateGenotypes(originalVC);
            if (regenotyped == null) {
                return null;
            }
            //make sure result has annotations so we don't have to keep originalVC around
            final Map<String, Object> newAnnotations;
            if (regenotyped.getNAlleles() != originalVC.getNAlleles()) {
                final Permutation<Allele> allelePermutation = new IndexedAlleleList<>(originalVC.getAlleles()).
                        permutation(new IndexedAlleleList<>(regenotyped.getAlleles()));
                final int[] relevantIndices = IntStream.range(0, regenotyped.getAlleles().size())
                        .map(n -> allelePermutation.fromIndex(n)).toArray();
                newAnnotations = new LinkedHashMap<>();
                composeUpdatedAnnotations(originalVC, newAnnotations, relevantIndices, regenotyped);
            } else {
                newAnnotations = originalVC.getAttributes();
            }
            result = new VariantContextBuilder(regenotyped).attributes(newAnnotations).make();
        }


        //variants with PL[0] less than threshold get turned to homRef with PL=[0,0,0], shouldn't get INFO attributes
        //make sure we can call het variants with GQ >= rgqThreshold in joint calling downstream
        if(shouldBeReblocked(result)) {
            if (result.getEnd() <= vcfOutputEnd) {
                return null;
            }
            final VariantContextBuilder newHomRefBuilder = lowQualVariantToGQ0HomRef(result);
            if (newHomRefBuilder != null) {
                updateHomRefBlockBuffer(newHomRefBuilder.make());
            }
            return null;  //don't write yet in case new ref block needs to be modified
        }
        //high quality variant
        else {
            final VariantContext trimmedVariant = cleanUpHighQualityVariant(result);

            //Handle overlapping deletions so there are no duplicate bases
            //Queue the low quality deletions until we hit a high quality variant or the start is past the oldBlockEnd
            updateHomRefBlockBuffer(trimmedVariant);
            return trimmedVariant;
        }
    }

    /**
     * Write and remove ref blocks that end before the variant
     * Trim ref block if the variant occurs in the middle of a block
     * @param variantContextToOutput can overlap existing ref blocks in buffer, but should never start before vcfOutputEnd
     */
    private void updateHomRefBlockBuffer(final VariantContext variantContextToOutput) {
        if (variantContextToOutput == null) {
            return;
        }
        Utils.validate(variantContextToOutput.getStart() <= variantContextToOutput.getEnd(),
                "Input variant context at position " + currentContig + ":" + variantContextToOutput.getStart() + " has negative length: start=" + variantContextToOutput.getStart() + " end=" + variantContextToOutput.getEnd());
        if (variantContextToOutput.getGenotype(0).isHomRef() && variantContextToOutput.getStart() < vcfOutputEnd) {
            throw new IllegalStateException("Reference positions added to buffer should not overlap positions already output to VCF. "
                    + variantContextToOutput.getStart() + " overlaps position " + currentContig + ":" + vcfOutputEnd + " already emitted.");
        }
        final List<VariantContextBuilder> completedBlocks = new ArrayList<>();
        final List<VariantContextBuilder> tailBuffer = new ArrayList<>();
        for (final VariantContextBuilder builder : homRefBlockBuffer) {
            final int blockStart = (int)builder.getStart();
            final int variantEnd = variantContextToOutput.getEnd();
            if (blockStart > variantEnd) {
                break;
            }
            int blockEnd = (int)builder.getStop();
            final int variantStart = variantContextToOutput.getStart();
            if (blockEnd >= variantStart) {  //then trim out overlap
                if (blockEnd > variantEnd && blockStart <= variantStart) {  //then this block will be split -- add a post-variant block
                    final VariantContextBuilder blockTailBuilder = new VariantContextBuilder(builder);
                    moveBuilderStart(blockTailBuilder, variantEnd + 1);
                    tailBuffer.add(blockTailBuilder);
                    builder.stop(variantEnd);
                    builder.attribute(VCFConstants.END_KEY, variantEnd);
                    blockEnd = variantEnd;
                }
                if (blockStart < variantStart) { //right trim
                    if (blockStart > variantStart - 1) {
                        throw new GATKException.ShouldNeverReachHereException("ref blocks screwed up; current builder: " + builder.getStart() + " to " + builder.getStop());
                    }
                    builder.attribute(VCFConstants.END_KEY, variantStart - 1);
                    builder.stop(variantStart - 1);
                } else {  //left trim
                    if (variantContextToOutput.contains(new SimpleInterval(currentContig, blockStart, blockEnd))) {
                        completedBlocks.add(builder);
                    } else {
                        if (blockEnd < variantEnd + 1) {
                            throw new GATKException.ShouldNeverReachHereException("ref blocks screwed up; current builder: " + builder.getStart() + " to " + builder.getStop());
                        }
                        moveBuilderStart(builder, variantEnd + 1);
                    }
                }
                if (builder.getStart() > builder.getStop()) {
                    throw new GATKException.ShouldNeverReachHereException("ref blocks screwed up; current builder: " + builder.getStart() + " to " + builder.getStop());
                }
            }
            //only flush ref blocks if we're outputting a variant, otherwise ref blocks can be out of order
            if (builder.getStart() < variantStart && !variantContextToOutput.getGenotype(0).isHomRef()) {
                vcfWriter.add(builder.make());
                vcfOutputEnd = (int)builder.getStop();
                completedBlocks.add(builder);
            }
            bufferEnd = blockEnd;  //keep track of observed ends
        }
        homRefBlockBuffer.removeAll(completedBlocks);
        final Genotype g = variantContextToOutput.getGenotype(0);
        if (g.isHomRef() || (g.hasPL() && g.getPL()[0] == 0)) {
            final VariantContextBuilder newHomRefBlock = new VariantContextBuilder(variantContextToOutput);
            homRefBlockBuffer.add(newHomRefBlock);
        }
        homRefBlockBuffer.addAll(tailBuffer);
        homRefBlockBuffer.sort(VCB_COMPARATOR);  //this may seem lazy, but it's more robust to assumptions about overlap being violated
        bufferEnd = Math.max(bufferEnd, variantContextToOutput.getEnd());
    }

    /**
     * Write all the reference blocks in the buffer to the output VCF
     */
    private void flushRefBlockBuffer() {
         for (final VariantContextBuilder builder : homRefBlockBuffer) {
             vcfWriter.add(builder.make());
             vcfOutputEnd = (int)builder.getStop();
         }
         homRefBlockBuffer.clear();
         bufferEnd = 0;
    }

    /**
     * determine if a VC is a homRef block, i.e. has an end key and does not have filtering annotations
     * @param result VariantContext to process
     * @return true if VC is a homRef block and not a "call" with annotations
     */
    private boolean isHomRefBlock(final VariantContext result) {
        if (result.getGenotype(0).hasPL()) {
            if (result.getGenotype(0).getPL()[0] != 0) {
                return false;
            }
        } else {
            return (result.getAttributes().size() == 1) && result.hasAttribute(VCFConstants.END_KEY);
        }
        return result.getLog10PError() == VariantContext.NO_LOG10_PERROR;
    }

    /**
     * determine if VC is a homRef "call", i.e. an annotated variant with non-symbolic alt alleles and homRef genotypes
     * we treat these differently from het/homVar calls or homRef blocks
     * @param result VariantContext to process
     * @return true if VC is a 0/0 call and not a homRef block
     */
    private boolean isMonomorphicCallWithAlts(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        return (hasGenotypeValuesArray(genotype)
                && (genotype.isHomRef() || isNoCalledHomRef(genotype) || (genotype.hasPL() && MathUtils.minElementIndex(genotype.getPL()) == 0))
                && result.getAlternateAlleles().stream().anyMatch(this::isConcreteAlt));
    }

    /**
     *
     * @param genotype  a genotype that may or may not have likelihood/probability data
     * @return true if we can quantify genotype call
     */
    private boolean hasGenotypeValuesArray(final Genotype genotype) {
        return genotype.hasPL() || (posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey));
    }

    private boolean isNoCalledHomRef(final Genotype genotype) {
        return genotype.isNoCall()
                && hasGenotypeValuesArray(genotype)
                && getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey)[0] == 0;
    }

    /**
     * Process a reference block VC:
     * Drop if low quality or covered by previous VCF output,
     * correct missing data if possible, move ref block start if necessary
     * @param result    a reference block
     * @return  null if this block is entirely overlapped by positions already written out or
     * if the ref block is low quality and we'redropping low quality variants
     */
    private VariantContext processHomRefBlock(final VariantContext result) {
        final Genotype genotype = result.getGenotype(0);
        final VariantContextBuilder vcBuilder = new VariantContextBuilder(result);
        if (result.getStart() <= vcfOutputEnd) {
            if (result.getEnd() <= vcfOutputEnd) {
                return null;
            }
            moveBuilderStart(vcBuilder, vcfOutputEnd + 1);
        }
        if (dropLowQuals && (!genotype.hasGQ() || genotype.getGQ() < rgqThreshold || genotype.getGQ() == 0)) {
            return null;
        } else if (genotype.isCalled() && genotype.isHomRef()) {
            if (!genotype.hasPL()) {
                if (genotype.hasGQ()) {
                    logger.warn("PL is missing for hom ref genotype at at least one position for sample " + genotype.getSampleName() + ": " + result.getContig() + ":" + result.getStart() +
                            ".  Using GQ to determine quality.");
                    final int gq = genotype.getGQ();
                    final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                    vcBuilder.genotypes(gBuilder.GQ(gq).make());
                    return vcBuilder.make();
                } else {
                    final String message = "Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at "
                            + result.getContig() + ":" + result.getStart() + " for sample " + genotype.getSampleName() + ".";
                    if (allowMissingHomRefData) {
                        logger.warn(message);
                        final GenotypeBuilder gBuilder = new GenotypeBuilder(genotype);
                        vcBuilder.genotypes(gBuilder.GQ(0).PL(new int[]{0,0,0}).make());
                        return vcBuilder.make();
                    } else {
                        throw new UserException.BadInput(message);
                    }
                }
            }
            return vcBuilder.make();
        //some external data has no-called genotypes with good likelihoods
        } else if (!genotype.isCalled() && genotype.hasPL() && genotype.getPL()[0] == 0) {
            return vcBuilder.make();
        }
        else {
            return null;
        }
    }

    /**
     * Should this variant context be turned into a reference block?
     * @param vc    a low quality variant or variant called homozygous reference
     * @return  true if this VariantContext is eligible to be combined with adjacent reference blocks
     */
    @VisibleForTesting
    boolean shouldBeReblocked(final VariantContext vc) {
        if (!vc.hasGenotypes()) {
            throw new IllegalStateException("Variant contexts must contain genotypes to be reblocked.");
        }
        final Genotype genotype = vc.getGenotype(0);
        final int[] pls = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
        return (pls != null && pls[0] < rgqThreshold) || genotype.isHomRef()
                || !genotypeHasConcreteAlt(genotype)
                || genotype.getAlleles().stream().anyMatch(a -> a.equals(Allele.NON_REF_ALLELE))
                || (!genotype.hasPL() && !genotype.hasGQ())
                || (genotype.hasDP() && genotype.getDP() == 0);
    }

    /**
     * Is there a real ALT allele that's not <NON_REF> or * ?
     * @param g called genotype
     * @return true if the genotype has a called allele that is a "real" alternate
     */
    private boolean genotypeHasConcreteAlt(final Genotype g) {
        return g.getAlleles().stream().anyMatch(this::isConcreteAlt);
    }

    private boolean isConcreteAlt(final Allele a) {
        return !a.isReference() && !a.isSymbolic() && !a.equals(Allele.SPAN_DEL);
    }

    /**
     * "reblock" a variant by converting its genotype to homRef, changing PLs, adding reblock END tags and other attributes
     * @param lowQualityVariant  a variant already determined to be low quality
     * @return a Builder that can be modified later
     */
    @VisibleForTesting
    public VariantContextBuilder lowQualVariantToGQ0HomRef(final VariantContext lowQualityVariant) {
        if(dropLowQuals && (!isMonomorphicCallWithAlts(lowQualityVariant) || !lowQualityVariant.getGenotype(0).isCalled())) {
            return null;
        }

        //change no-calls with no PL values to ref blocks
        if (genotypeHasNoData(lowQualityVariant.getGenotype(0))) {
            final Map<String, Object> blockAttributes = new LinkedHashMap<>(2);
            final GenotypeBuilder gb = changeCallToHomRefVersusNonRef(lowQualityVariant, blockAttributes);
            final List<Allele> blockAlleles = Arrays.asList(Allele.create(lowQualityVariant.getReference().getBases()[0], true), Allele.NON_REF_ALLELE);
            final VariantContextBuilder vb = new VariantContextBuilder(lowQualityVariant).alleles(blockAlleles).attributes(blockAttributes).genotypes(gb.make())
                    .log10PError(VariantContext.NO_LOG10_PERROR);
            if (vcfOutputEnd + 1 > lowQualityVariant.getStart()) {
                moveBuilderStart(vb, vcfOutputEnd + 1);
            }
            return vb;   //delete variant annotations and add end key
        }

        final Map<String, Object> attrMap = new HashMap<>();

        //this method does a lot of things, including fixing alleles and adding the END key
        final GenotypeBuilder gb = changeCallToHomRefVersusNonRef(lowQualityVariant, attrMap);  //note that gb has all zero PLs

        final VariantContextBuilder builder = new VariantContextBuilder(lowQualityVariant);

        final Genotype newG = gb.make();
        builder.alleles(Arrays.asList(newG.getAlleles().get(0), Allele.NON_REF_ALLELE)).genotypes(newG);
        if (lowQualityVariant.getStart() <= vcfOutputEnd) {
            final int newStart = vcfOutputEnd + 1;
            moveBuilderStart(builder, newStart);
        }
        return builder.unfiltered()
                .log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap); //genotyping engine will add lowQual filter, so strip it off
    }

    private boolean genotypeHasNoData(final Genotype g) {
        return !g.isCalled() && !g.hasPL();
    }

    /**
     * Produce a GenotypeBuilder for a hom-ref call suitable to be merged into a reference block, i.e. set with PLs as
     * appropriate for a variant with only reference and <NON_REF> as alleles and END applied as appropriate
     * Note that this may modify {@code attrMap} as a side effect for END key
     * @param lowQualVariant a VC containing a genotype to be converted to a GQ0 homRef call; needed for alleles that correspond to PLs and other things
     * @param attrMap the new VC attribute map, to update the END tag for deletions
     * @return a GenotypeBuilder to make a 0/0 call with updated PLs and GQ
     */
    @VisibleForTesting
    protected GenotypeBuilder changeCallToHomRefVersusNonRef(final VariantContext lowQualVariant, final Map<String, Object> attrMap) {
        final Genotype genotype = lowQualVariant.getGenotype(0);
        final Allele inputRefAllele = lowQualVariant.getReference();
        GenotypeBuilder gb = new GenotypeBuilder(genotype);
        //if GT is not homRef, correct it and set GQ=0
        if (posteriorsKey == null && (!genotype.hasPL() || (genotype.hasPL() && genotype.getPL()[0] != 0))) {
            gb.PL(new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())]);  //2 alleles for ref and non-ref
            gb.GQ(0).noAD().alleles(Collections.nCopies(genotype.getPloidy(), inputRefAllele)).noAttributes();
        //for hom-ref variants, drop other ALTs and subset PLs, GQ is recalculated (may be > 0)
        } else {
            if (posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey)) {
                subsetPosteriorsToRefVersusNonRef(lowQualVariant, gb);
            } else {
                final List<Allele> bestAlleles = AlleleSubsettingUtils.calculateMostLikelyAllelesForMonomorphicSite(lowQualVariant, PLOIDY_TWO, 1);
                final Allele bestAlt = bestAlleles.stream().filter(a -> !a.isReference() && !a.isNonRefAllele()).findFirst().orElse(Allele.NON_REF_ALLELE);  //allow span dels
                //here we're assuming that an alt that isn't <NON_REF> will have higher likelihood than non-ref, which should be true
                final GenotypesContext context = AlleleSubsettingUtils.subsetAlleles(lowQualVariant.getGenotypes(),
                        genotype.getPloidy(), lowQualVariant.getAlleles(), Arrays.asList(inputRefAllele, bestAlt),
                        null, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, lowQualVariant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));  //BEST_MATCH to avoid no-calling low qual genotypes
                final Genotype subsetG = context.get(0);
                gb = new GenotypeBuilder(subsetG).noAttributes();  //remove attributes because hom ref blocks shouldn't have posteriors
                //subsetting may strip GQ and PLs for low qual genotypes
                if (!subsetG.hasGQ()) {
                    gb.GQ(0);
                }
                if (!subsetG.hasPL()) {
                    gb.PL(new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())]);  //2 alleles for ref and non-ref
                }
            }
        }
        if (lowQualVariant.hasAttribute(VCFConstants.DEPTH_KEY)) {
            final int depth = lowQualVariant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            gb.DP(depth);
            gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, depth);
        } else if (genotype.hasAD()) {
            final int depth = (int) MathUtils.sum(genotype.getAD());
            gb.DP(depth);
            gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, depth);
        }
        //NOTE: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
        final Allele outputRefAllele;
        if (inputRefAllele.length() > 1 || genotype.getAlleles().contains(Allele.SPAN_DEL) || genotype.getAlleles().contains(Allele.NO_CALL)) {
            outputRefAllele = Allele.create(inputRefAllele.getBases()[0], true);
        } else {
            outputRefAllele = inputRefAllele;
        }
        attrMap.put(VCFConstants.END_KEY, lowQualVariant.getEnd());
        gb.alleles(Collections.nCopies(genotype.getPloidy(), outputRefAllele));
        return gb;
    }

    /**
     * Subset alleles as necessary and apply annotations
     * @param variant    VariantContext with full set of annotations (e.g. DP)
     * @return  an annotated VariantContext with data only for ref, non-ref and called alts
     */
    @VisibleForTesting
    VariantContext cleanUpHighQualityVariant(final VariantContext variant) {
        final Map<String, Object> attrMap = new HashMap<>();

        final Genotype genotype = getCalledGenotype(variant);
        VariantContextBuilder builder = new VariantContextBuilder(variant);  //QUAL from result is carried through
        builder.attributes(attrMap);  //clear attributes

        final List<Allele> allelesToDrop = getAllelesToDrop(variant, genotype);

        final boolean allelesNeedSubsetting = !allelesToDrop.isEmpty();
        int[] relevantIndices = new int[variant.getNAlleles()];  //called alleles plus ref and non-ref
        final List<Allele> newAlleleSetUntrimmed = new ArrayList<>(variant.getAlleles());
        if(allelesNeedSubsetting && !keepAllAlts) {
            newAlleleSetUntrimmed.removeAll(allelesToDrop);
            final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(variant.getGenotypes(), PLOIDY_TWO, variant.getAlleles(),
                    newAlleleSetUntrimmed, null, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN,
                    variant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            if (!gc.get(0).hasGQ()) {
                logger.warn("No GQ returned for genotype at " + variant.getContig() + ":"
                        + variant.getStart() + " after subsetting -- is this really a high quality genotype? " + gc.get(0).toString());
                final VariantContextBuilder newHomRefBuilder = lowQualVariantToGQ0HomRef(variant);
                if (newHomRefBuilder != null) {
                    updateHomRefBlockBuffer(newHomRefBuilder.make());
                }
            }
            //note that subsetting alleles can increase GQ, e.g. with one confident reference allele and a deletion allele that's either 4 or 5 bases long
            builder.genotypes(gc).alleles(newAlleleSetUntrimmed);
            //if deletions are dropped, alleles may need trimming
            final VariantContext newTrimmedAllelesVC = GATKVariantContextUtils.trimAlleles(builder.make(), false, true);
            builder = new VariantContextBuilder(newTrimmedAllelesVC);
            //save indices of new alleles for annotation processing
            relevantIndices = newAlleleSetUntrimmed.stream().mapToInt(a -> variant.getAlleles().indexOf(a)).toArray();
            final int refBlockDepth;
            if (variant.hasAttribute(VCFConstants.DEPTH_KEY)) {  //prefer INFO depth because HaplotypeCaller GVCF block code uses all reads, not just informative
                refBlockDepth = variant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            } else if (genotype.hasDP()) {
                refBlockDepth = genotype.getDP();
            } else {
                refBlockDepth = 0;
            }
            addRefBlockIfNecessary(variant, allelesToDrop, newTrimmedAllelesVC, refBlockDepth);
        }

        final VariantContext updatedAllelesVC = builder.make();
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);

        //remove any AD reads for the non-ref
        final ArrayList<Genotype> genotypesArray = removeNonRefADs(updatedAllelesGenotype, updatedAllelesVC.getAlleleIndex(Allele.NON_REF_ALLELE));
        builder.genotypes(genotypesArray);

        composeUpdatedAnnotations(variant, attrMap, relevantIndices, updatedAllelesVC);

        return builder.attributes(attrMap).unfiltered().make();
    }

    /**
     * Update annotations if alleles were subset and add annotations specific to ReblockGVCF, like QUALapprox and RAW_GT_COUNT
     * @param variant   variant context with full annotation data
     * @param attrMap   annotation map to modify with new annotations
     * @param relevantIndices   indexes of alleles in updatedAllelesVC with respect to variant
     * @param updatedAllelesVC  variant context with final set of alleles
     */
    private void composeUpdatedAnnotations(final VariantContext variant, final Map<String, Object> attrMap,
                                           final int[] relevantIndices, final VariantContext updatedAllelesVC) {
        updateMQAnnotations(variant, attrMap);

        final boolean allelesNeedSubsetting = relevantIndices.length < variant.getNAlleles();
        copyInfoAnnotations(variant, attrMap, allelesNeedSubsetting, relevantIndices);

        //generate qual annotations after we potentially drop alleles
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);
        if (doQualApprox) {
            //TODO: update to support DRAGEN posteriors
            if (updatedAllelesGenotype.hasPL()) {
                addQualAnnotations(attrMap, updatedAllelesVC);
            }
        } else {  //manually copy annotations that might be from reblocking and aren't part of AnnotationEngine
            if (variant.hasAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY)) {
                attrMap.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, variant.getAttribute(GATKVCFConstants.AS_VARIANT_DEPTH_KEY));
            }
            if (variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
                attrMap.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, variant.getAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
            }
        }
        attrMap.put(GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY, updatedAllelesGenotype.getAlleles().stream().anyMatch(Allele::isReference) ?
                Arrays.asList(0,1,0) : Arrays.asList(0,0,1)); //ExcessHet currently uses rounded/integer genotype counts, so do the same here
    }

    /**
     * Return the genotype from variant, call it if data is present and GT is no-call
     * @param variant   VariantContext with genotype data
     * @return  a called genotype (if possible)
     */
    private Genotype getCalledGenotype(final VariantContext variant) {
        if (variant.getGenotype(0).isNoCall()) {
            final Genotype noCallGT = variant.getGenotype(0);
            final GenotypeBuilder builderToCallAlleles = new GenotypeBuilder(noCallGT);
            //TODO: update to support DRAGEN posteriors
            GATKVariantContextUtils.makeGenotypeCall(noCallGT.getPloidy(), builderToCallAlleles, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN,
                    noCallGT.getLikelihoods().getAsVector(), variant.getAlleles(), null);
            return builderToCallAlleles.make();
        } else {
            return variant.getGenotype(0);
        }
    }

    /**
     * Get the list of concrete alternate alleles in variant that were not called in genotype
     * @param variant   has full set of alleles, not all of which may be called
     * @param calledGenotype    should not be no-call
     * @return  a list of (concrete) alt alleles that are not called in calledGenotype
     */
    private List<Allele> getAllelesToDrop(final VariantContext variant, final Genotype calledGenotype) {
        //always drop alleles that aren't called to reduce PL size
        final List<Allele> allelesToDrop = new ArrayList<>();
        for (final Allele currAlt : variant.getAlternateAlleles()) {
            boolean foundMatch = false;
            for (final Allele gtAllele : calledGenotype.getAlleles()) {
                if (gtAllele.equals(currAlt, false)) {
                    foundMatch = true;
                    break;
                }
            }
            if (!foundMatch && !currAlt.isSymbolic()) {
                allelesToDrop.add(currAlt);
            }
        }
        return allelesToDrop;
    }

    /**
     * Add the "raw" annotations necessary for calculating QD and AS_QD
     * @param attrMap   has qual-related annotations added to it, but also potentially supplies DP value
     * @param updatedAllelesVC  variant context without uncalled alts
     */
    private void addQualAnnotations(final Map<String, Object> attrMap, final VariantContext updatedAllelesVC) {
        final Genotype updatedAllelesGenotype = updatedAllelesVC.getGenotype(0);
        attrMap.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, updatedAllelesGenotype.getPL()[0]);
        int varDP = QualByDepth.getDepth(updatedAllelesVC.getGenotypes(), null);
        if (varDP == 0) {  //prevent QD=Infinity case
            //attrMap should have DP already from copyInfoAnnotations call above
            varDP = Integer.parseInt(attrMap.getOrDefault(VCFConstants.DEPTH_KEY, 1).toString()); //if there's no VarDP and no DP, just prevent Infs/NaNs and QD will get capped later
        }
        attrMap.put(GATKVCFConstants.VARIANT_DEPTH_KEY, varDP);
        if (annotationEngine.getInfoAnnotations().stream()
                .anyMatch(infoFieldAnnotation -> infoFieldAnnotation.getClass().getSimpleName().equals("AS_QualByDepth"))) {
            final List<String> quals = new ArrayList<>();
            //get allele-specific QUAL approximation by subsetting PLs for each alt
            for (final Allele alt : updatedAllelesVC.getAlleles()) {
                if (alt.isReference()) {
                    //GDB expects an empty field for ref
                    continue;
                }
                if (alt.equals(Allele.NON_REF_ALLELE) || alt.equals(Allele.SPAN_DEL)) {
                    quals.add("0");
                    continue;
                }
                //TODO: this isn't going to work for DRAGEN's genotype posteriors
                final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(updatedAllelesVC.getGenotypes(),
                        updatedAllelesGenotype.getPloidy(), updatedAllelesVC.getAlleles(), Arrays.asList(updatedAllelesVC.getReference(), alt), null,
                        GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, 0);  //don't need depth to get PLs for quals

                final Genotype subsettedGenotype = gc.get(0);
                final int[] likelihoods = getGenotypePosteriorsOtherwiseLikelihoods(subsettedGenotype, posteriorsKey);
                if (likelihoods != null) {
                    quals.add(Integer.toString(likelihoods[0]));
                }  else {  //AlleleSubsettingUtils can no-call genotypes with super duper low GQs
                    quals.add("0");
                }
            }
            attrMap.put(GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY, AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM + String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, quals));
            final List<Integer> as_varDP = AS_QualByDepth.getAlleleDepths(updatedAllelesVC.getGenotypes());
            if (as_varDP != null) {
                attrMap.put(GATKVCFConstants.AS_VARIANT_DEPTH_KEY, as_varDP.stream().map(n -> Integer.toString(n)).collect(Collectors.joining(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM)));
            }
        }
    }

    /**
     * Any AD counts for the non-ref allele get propagated to every new allele when GVCFs are merged, so zero them out
     * @param g a genotype that may or may not contain AD
     * @param nonRefInd allele index of the non-ref, -1 if missing
     * @return  a Genotype array that can be used by a GenotypeBuilder
     */
    private ArrayList<Genotype> removeNonRefADs(final Genotype g, final int nonRefInd) {
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        if (g.hasAD() && nonRefInd != -1) {
            final int[] ad = g.getAD();
            if (ad.length >= nonRefInd && ad[nonRefInd] > 0) { //only initialize a builder if we have to
                final GenotypeBuilder gb = new GenotypeBuilder(g);
                ad[nonRefInd] = 0;
                gb.AD(ad).DP((int) MathUtils.sum(ad));
                genotypesArray.add(gb.make());
            } else {
                genotypesArray.add(g);
            }
        } else {
            genotypesArray.add(g);
        }
        return genotypesArray;
    }

    /**
     * Add the original annotations to the map for the new VC, subsetting AS annotations as necessary
     * @param originalVC    VC with full set of alleles and INFO annotations
     * @param newVCAttrMap   map of new annotations, to be modified
     * @param allelesNeedSubsetting do we have to subset allele-specific annotations?
     * @param relevantIndices   indexes for called alleles within the full alleles set from original VC
     */
    private void copyInfoAnnotations(final VariantContext originalVC, final Map<String, Object> newVCAttrMap,
                                     final boolean allelesNeedSubsetting, final int[] relevantIndices) {
        //copy over info annotations
        final Map<String, Object> origMap = originalVC.getAttributes();
        for(final InfoFieldAnnotation annotation : annotationEngine.getInfoAnnotations()) {
            for (final String key : annotation.getKeyNames()) {
                if (infoFieldAnnotationKeyNamesToRemove.contains(key)) {
                    continue;
                }
                if (origMap.containsKey(key)) {
                    newVCAttrMap.put(key, origMap.get(key));
                }
            }
            if (annotation instanceof ReducibleAnnotation) {
                for (final String rawKey : ((ReducibleAnnotation)annotation).getRawKeyNames()) {
                    if (infoFieldAnnotationKeyNamesToRemove.contains(rawKey)) {
                        continue;
                    }
                    if (origMap.containsKey(rawKey)) {
                        if (allelesNeedSubsetting && AnnotationUtils.isAlleleSpecific(annotation)) {
                            final List<String> alleleSpecificValues = AnnotationUtils.getAlleleLengthListOfString(originalVC.getAttributeAsString(rawKey, null));
                            final List<?> subsetList = alleleSpecificValues.size() > 0 ? AlleleSubsettingUtils.remapRLengthList(alleleSpecificValues, relevantIndices, "")
                                    : Collections.nCopies(relevantIndices.length, "");
                            newVCAttrMap.put(rawKey, AnnotationUtils.encodeAnyASListWithRawDelim(subsetList));
                        } else {
                            newVCAttrMap.put(rawKey, origMap.get(rawKey));
                        }
                    }
                }
            }
        }
    }

    /**
     * Add the newest raw mapping quality annotations to the annotation map
     * @param originalVC    VariantContext that may contain deprecated mapping quality annotations
     * @param newVCAttrMap   modified to add mapping quality raw annotations
     */
    private void updateMQAnnotations(final VariantContext originalVC, final Map<String, Object> newVCAttrMap) {
        //all VCs should get new RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, but preserve deprecated keys if present
        if (!originalVC.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY)) {
            //we're going to approximate depth for MQ calculation with the site-level DP (should be informative and uninformative reads),
            //which is pretty safe because it will only differ if reads are missing MQ
            final Integer rawMqValue = originalVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED) ?
                    (int)Math.round(originalVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0.0)) :
                    (int)Math.round(originalVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            originalVC.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, 60.0) *
                            originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            newVCAttrMap.put(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY,
                    String.format("%d,%d", rawMqValue, originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0)));

            //NOTE: this annotation is deprecated, but keep it here so we don't have to reprocess older GVCFs
            if (originalVC.hasAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
                newVCAttrMap.put(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED,
                        originalVC.getAttributeAsDouble(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0));
                newVCAttrMap.put(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED,
                        originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            }
        }
    }

    /**
     * If the ref allele is trimmed after alt deletions are dropped, add a reference block to account for the space covered before trimming
     * @param result    VC with full set of alleles that may need to be trimmed
     * @param allelesToDrop alleles eligible to become the new non-ref likelihood
     * @param newTrimmedAllelesVC   VC with called alleles that may have been trimmed
     */
    private void addRefBlockIfNecessary(final VariantContext result, final List<Allele> allelesToDrop, final VariantContext newTrimmedAllelesVC, final int refBlockDepth) {
        //if deletion needs trimming, fill in the gap with a ref block
        final int oldLongestAlleleLength = result.getReference().length();
        final int newLongestAlleleLength = newTrimmedAllelesVC.getReference().length();
        final Genotype genotype = result.getGenotype(0);
        if (newLongestAlleleLength < oldLongestAlleleLength) {
            //need to add a ref block to make up for the allele trimming or there will be a hole in the GVCF
            final int[] originalLikelihoods = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
            if (originalLikelihoods != null) {
                final Allele oldShortestAltAllele;
                try {
                    oldShortestAltAllele = allelesToDrop.stream().filter(a -> !a.equals(Allele.SPAN_DEL)).min(Allele::compareTo).orElseThrow(NoSuchElementException::new);
                } catch (final Exception e) {
                    throw new GATKException("No shortest ALT at " + result.getStart() + " across alleles: " + allelesToDrop);
                }

                //subset PLs to ref and longest dropped allele (longest may not be most likely, but we'll approximate so we don't have to make more than one ref block)
                final int[] longestVersusRefPLIndices = AlleleSubsettingUtils.subsettedPLIndices(result.getGenotype(0).getPloidy(),
                        result.getAlleles(), Arrays.asList(result.getReference(), oldShortestAltAllele));
                final int[] newRefBlockLikelihoods = MathUtils.normalizePLs(Arrays.stream(longestVersusRefPLIndices)
                        .map(idx -> originalLikelihoods[idx]).toArray());
                if (newRefBlockLikelihoods[0] != 0) {
                    for (int i = 0; i < newRefBlockLikelihoods.length; i++) {
                        newRefBlockLikelihoods[i] = Math.max(newRefBlockLikelihoods[i] - newRefBlockLikelihoods[0], 0);
                    }
                }

                //build the new reference block with updated likelihoods
                final GenotypeBuilder refBlockGenotypeBuilder = new GenotypeBuilder();
                final int refStart = Math.max(result.getEnd() - (oldLongestAlleleLength - newLongestAlleleLength), vcfOutputEnd) + 1;
                final Allele newRef = Allele.create(ReferenceUtils.getRefBaseAtPosition(referenceReader, result.getContig(), refStart), true);
                refBlockGenotypeBuilder.PL(newRefBlockLikelihoods)
                        .GQ(MathUtils.secondSmallestMinusSmallest(newRefBlockLikelihoods, 0))
                        .alleles(Arrays.asList(newRef, newRef)).DP(refBlockDepth);

                //add the new block to the buffer if it isn't covered by positions already output
                if (refStart > vcfOutputEnd && result.getEnd() > vcfOutputEnd) {
                    final VariantContextBuilder trimBlockBuilder = new VariantContextBuilder();
                    trimBlockBuilder.chr(currentContig).start(Math.max(refStart, vcfOutputEnd + 1)).stop(result.getEnd()).
                            alleles(Arrays.asList(newRef, Allele.NON_REF_ALLELE)).attribute(VCFConstants.END_KEY, result.getEnd())
                            .genotypes(refBlockGenotypeBuilder.make());
                    updateHomRefBlockBuffer(trimBlockBuilder.make());
                }
            }
        }
    }

    /**
     * Modifies ref block builder to change start position and update ref allele accordingly in VC and genotypes
     * @param builder   a builder for a reference block
     * @param newStart  the new position for the reference block
     */
    private void moveBuilderStart(final VariantContextBuilder builder, final int newStart) {
        final byte[] newRef = ReferenceUtils.getRefBaseAtPosition(referenceReader, currentContig, newStart);
        final Allele newRefAllele = Allele.create(newRef, true);
        final ArrayList<Genotype> genotypesArray = new ArrayList<>();
        for (final Genotype g : builder.getGenotypes()) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);
            final List<Allele> newGTAlleles = g.getAlleles().stream().map(a -> a.isReference() ? newRefAllele : a).collect(Collectors.toList());
            gb.alleles(newGTAlleles);
            genotypesArray.add(gb.make());
        }
        final List<Allele> newVCAlleles = builder.getAlleles().stream().map(a -> a.isReference() ? newRefAllele : a).collect(Collectors.toList());
        builder.start(newStart).alleles(newVCAlleles).genotypes(genotypesArray);
    }

    /**
     * If genotype posterior probabilities are present, return those; otherwise use likelihoods
     * @param genotype  should have PLs
     * @param posteriorsKey may be null
     * @return may be null
     */
    private int[] getGenotypePosteriorsOtherwiseLikelihoods(final Genotype genotype, final String posteriorsKey) {
        if ((posteriorsKey != null && genotype.hasExtendedAttribute(posteriorsKey))) {
            final double[] posteriors = VariantContextGetters.getAttributeAsDoubleArray(genotype, posteriorsKey, () -> null, 0);
            return Arrays.stream(posteriors).mapToInt(x -> (int)Math.round(x)).toArray();
        } else if (genotype.hasPL()) {
            return genotype.getPL();
        } else {
            return null;
        }
    }

    /**
     * Given a variant with multi-allelic PLs, modify the GenotypeBuilder to have annotations as for just ref and non-ref
     * @param result    a VariantContext containing alternate alleles in addition to non-ref
     * @param gb    a reference block GenotypeBuilder to be modified
     */
    private void subsetPosteriorsToRefVersusNonRef(final VariantContext result, final GenotypeBuilder gb) {
        //TODO: bestAlleles needs to be modified for posteriors
        final List<Allele> bestAlleles = AlleleSubsettingUtils.calculateMostLikelyAlleles(result, PLOIDY_TWO, 1);
        final Allele bestAlt = bestAlleles.stream().filter(a -> !a.isReference() && !a.isNonRefAllele()).findFirst().orElse(Allele.NON_REF_ALLELE);  //allow span dels
        //here we're assuming that an alt that isn't <NON_REF> will have higher likelihood than non-ref, which should be true
        final int[] idxVector = result.getGLIndicesOfAlternateAllele(bestAlt);
        final Genotype genotype = result.getGenotype(0);
        final int[] multiallelicPLs = getGenotypePosteriorsOtherwiseLikelihoods(genotype, posteriorsKey);
        if (multiallelicPLs != null) {
            int[] newPLs = new int[GenotypeLikelihoods.numLikelihoods(2, genotype.getPloidy())];
            for (int i = 0; i < idxVector.length; i++) {
                newPLs[i] = multiallelicPLs[idxVector[i]];
            }
            //in the case of *, we need to renormalize to homref
            if (newPLs[0] != 0) {
                final int[] output = new int[newPLs.length];
                for (int i = 0; i < newPLs.length; i++) {
                    output[i] = Math.max(newPLs[i] - newPLs[0], 0);
                }
                newPLs = output;
            }
            gb.PL(newPLs);
            gb.GQ(MathUtils.secondSmallestMinusSmallest(newPLs, 0));
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
