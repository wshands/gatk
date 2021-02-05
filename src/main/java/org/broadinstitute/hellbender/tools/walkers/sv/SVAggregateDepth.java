package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.sv.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Calculates copy number posteriors for a given set of structural variants.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 *     <li>
 *         Germline CNV intervals VCF
 *     </li>
 *     <li>
 *         Germline CNV contig ploidy calls
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF with copy number posteriors
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVAggregateDepth
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Collects read counts at specified intervals",
        oneLineSummary = "Collects read counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class SVAggregateDepth extends VariantWalker {
    public static final String COPY_NUMBER_INTERVALS_LONG_NAME = "cnv-intervals-vcf";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls-file";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String COPY_NEUTRAL_PRIOR_LONG_NAME = "copy-neutral-prior";
    public static final String CNV_BND_PROB_THRESHOLD_LONG_NAME = "cnv-bnd-max-phred";
    public static final String CNV_BND_SAMPLE_THRESHOLD_LONG_NAME = "cnv-bnd-carrier-fraction";
    public static final String GENOTYPE_DEPTH_CALLS_LONG_NAME = "genotype-depth-calls";

    @Argument(
            doc = "Germline copy number intervals VCF. Can be specified more than once for runs with different bin sizes.",
            fullName = COPY_NUMBER_INTERVALS_LONG_NAME
    )
    private List<File> copyNumberPosteriorsFiles;

    @Argument(
            doc = "Contig ploidy calls file. Can be specified for multiple samples.",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private List<File> contigPloidyCallFiles;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Perform genotyping on depth-only calls.",
            fullName = GENOTYPE_DEPTH_CALLS_LONG_NAME
    )
    private boolean genotypeDepthCalls = false;

    private List<VCFFileReader> posteriorsReaders;
    private List<String> samples;
    private VariantContextWriter outputWriter;
    private SVGenotypeEngineDepthOnly depthOnlyGenotypeEngine;
    private DepthEvidenceAggregator depthEvidenceAggregator;
    private SAMSequenceDictionary dictionary;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        if (genotypeDepthCalls) {
            depthOnlyGenotypeEngine = new SVGenotypeEngineDepthOnly();
        }
        posteriorsReaders = copyNumberPosteriorsFiles.stream().map(VCFFileReader::new).collect(Collectors.toList());
        samples = getHeaderForVariants().getSampleNamesInOrder();
        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFiles.stream()
                .map(CalledContigPloidyCollection::new)
                .collect(Collectors.toList());
        validateSampleSets(contigPloidyCollections);
        depthEvidenceAggregator = new DepthEvidenceAggregator(posteriorsReaders, contigPloidyCollections, samples, dictionary);
        outputWriter = createVCFWriter(outputFile);
        outputWriter.writeHeader(composeHeader());
    }

    @Override
    public Object onTraversalSuccess() {
        outputWriter.close();
        return null;
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecordUtils.create(variant);
        final SVCallRecordDepthPosterior posterior = depthEvidenceAggregator.apply(call);
        VariantContext finalVariant = posterior == null ? variant : createBaseVariant(variant, posterior);
        if (genotypeDepthCalls && SVGenotypeEngineFromModel.isDepthOnlyVariant(finalVariant)) {
            finalVariant = depthOnlyGenotypeEngine.genotypeVariant(finalVariant);
        }
        outputWriter.add(finalVariant);
    }

    private void validateSampleSets(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        final Set<String> svSamples = new HashSet<>(samples);
        final Set<String> cnvSamples = getCNVSamples();
        final Set<String> contigPloidySamples = getContigPloidySamples(contigPloidyCollections);
        validateSubset(cnvSamples, svSamples, "The following samples from the SV VCF were not in the CNV VCF");
        validateSubset(contigPloidySamples, svSamples, "The following samples from the SV VCF were not in the contig ploidy calls");
    }

    private void validateSubset(final Set<String> set, final Set<String> subset, final String msg) {
        final Set<String> diff = Sets.difference(subset, set);
        Utils.validate(diff.isEmpty(), msg + ":" + String.join(", ", diff));
    }

    private Set<String> getContigPloidySamples(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        return contigPloidyCollections.stream()
                .map(p -> p.getMetadata().getSampleName())
                .collect(Collectors.toSet());
    }

    private Set<String> getCNVSamples() {
        final Set<String> cnvSamplesSet = new HashSet<>(posteriorsReaders.get(0).getFileHeader().getSampleNamesInOrder());
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final Set<String> set = new HashSet<>(posteriorsReaders.get(i).getFileHeader().getSampleNamesInOrder());
            if (!set.equals(cnvSamplesSet)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical samples.");
            }
        }
        return Collections.unmodifiableSet(cnvSamplesSet);
    }

    private VCFHeader composeHeader() {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(getDefaultToolVCFHeaderLines());
        headerInfo.add(new VCFFormatHeaderLine(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, 1,
                VCFHeaderLineType.Integer, "Phred-scaled copy number posterior over the event region"));
        headerInfo.add(new VCFFormatHeaderLine(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY, 1,
                VCFHeaderLineType.Integer, "Neutral copy number"));
        headerInfo.add(new VCFInfoHeaderLine(GATKSVVCFConstants.DEPTH_OVERLAP_KEY, 1,
                VCFHeaderLineType.Integer, "Length of variant overlapping read depth evidence (CNV/BND only)"));
        if (genotypeDepthCalls) {
            headerInfo.addAll(depthOnlyGenotypeEngine.getHeaderLines());
        }
        final VCFHeader vcfHeader = new VCFHeader(getHeaderForVariants());
        headerInfo.stream().forEach(line -> vcfHeader.addMetaDataLine(line));
        return vcfHeader;
    }

    private VariantContext createBaseVariant(final VariantContext variant, final SVCallRecordDepthPosterior posteriors) {
        Utils.nonNull(variant);
        Utils.nonNull(posteriors);
        final VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        variantBuilder.attribute(GATKSVVCFConstants.DEPTH_OVERLAP_KEY, posteriors.getPosteriorIntervalsOverlap());
        final Map<String, List<Integer>> copyStateQuals = getCopyStateQuals(posteriors.getSamplePosteriors());
        final List<Genotype> newGenotypes = new ArrayList<>(variant.getGenotypes().size());
        for (final Genotype genotype : variant.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            final String sample = genotype.getSampleName();
            genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, copyStateQuals.get(sample));
            genotypeBuilder.attribute(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY, depthEvidenceAggregator.getSamplePloidy(sample, variant.getContig()));
            newGenotypes.add(genotypeBuilder.make());
        }
        variantBuilder.genotypes(newGenotypes);
        return variantBuilder.make();
    }

    private Map<String, List<Integer>> getCopyStateQuals(final Map<String, CopyNumberPosteriorDistribution> variantPosteriors) {
        final Map<String, List<Integer>> copyStateQuals = new HashMap<>(SVUtils.hashMapCapacity(variantPosteriors.size()));
        for (final String sample : variantPosteriors.keySet()) {
            final CopyNumberPosteriorDistribution dist = variantPosteriors.get(sample);
            final List<Integer> copyStatePhred = dist.getIntegerCopyNumberStateList().stream()
                    .map(dist::getCopyNumberPosterior)
                    .map(QualityUtils::logProbToPhred)
                    .map(Math::round)
                    .map(Long::intValue)
                    .collect(Collectors.toList());
            copyStateQuals.put(sample, copyStatePhred);
        }
        return copyStateQuals;
    }
}
