package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public class DepthEvidenceAggregator {

    private final List<VCFFileReader> posteriorsReaders;
    private final SAMSequenceDictionary dictionary;
    private final List<String> samples;
    private final List<IntegerCopyNumberState> copyStates;
    private final int numCopyStates;
    private Map<String,Map<String,Integer>> sampleContigPloidyMap;

    private String currentContig;
    private List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList;

    public DepthEvidenceAggregator(final List<VCFFileReader> posteriorsReaders,
                                   final Collection<CalledContigPloidyCollection> contigPloidyCollections,
                                   final List<String> samples,
                                   final SAMSequenceDictionary dictionary) {
        Utils.nonNull(posteriorsReaders);
        Utils.nonNull(contigPloidyCollections);
        Utils.nonNull(samples);
        Utils.nonNull(dictionary);
        this.posteriorsReaders = posteriorsReaders;
        this.samples = samples;
        this.dictionary = dictionary;
        this.sampleContigPloidyMap = getSampleContigPloidyMap(contigPloidyCollections);
        final VariantContext exampleVariant = posteriorsReaders.get(0).iterator().next();
        copyStates = getCopyNumberStates(exampleVariant);
        numCopyStates = copyStates.size();
        validateCopyStates();
    }

    public int getSamplePloidy(final String sample, final String contig) {
        if (!sampleContigPloidyMap.containsKey(sample)) {
            throw new UserException("Could not find ploidy calls for sample: " + sample);
        }
        final Map<String,Integer> contigPloidyMap = sampleContigPloidyMap.get(sample);
        if (!contigPloidyMap.containsKey(contig)) {
            throw new UserException("Could not find ploidy calls for contig: " + contig);
        }
        return contigPloidyMap.get(contig);
    }

    public SVCallRecordDepthPosterior apply(final SVCallRecord call) {
        Utils.nonNull(call);
        if (!call.getContigA().equals(currentContig)) {
            currentContig = call.getContigA();
            currentPosteriorsTreeList = posteriorsReaders.stream().map(this::getCurrentPosteriorsTree).collect(Collectors.toList());
        }
        if (!(call.getType().equals(StructuralVariantType.DEL)  || call.getType().equals(StructuralVariantType.DUP)
                || call.getType().equals(StructuralVariantType.BND))) {
            return null;
        }
        if (!call.getContigA().equals(call.getContigB())) {
            return null;
        }
        return SVCallRecordDepthPosterior.create(copyStates, samples, call, currentPosteriorsTreeList);
    }

    private void validateCopyStates() {
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final List<IntegerCopyNumberState> otherCopyStates = getCopyNumberStates(posteriorsReaders.get(i).iterator().next());
            if (!copyStates.equals(otherCopyStates)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical copy number states.");
            }
        }
    }

    private List<IntegerCopyNumberState> getCopyNumberStates(final VariantContext variant) {
        final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
        return parseGenotype(variant.getGenotype(0), interval)
                .getCopyNumberPosteriorDistribution().getIntegerCopyNumberStateList();
    }

    private IntervalTree<Map<String,double[]>> getCurrentPosteriorsTree(final VCFFileReader reader) {
        final SAMSequenceRecord contigRecord = dictionary.getSequence(currentContig);
        if (contigRecord == null) {
            throw new UserException.MissingContigInSequenceDictionary(currentContig, dictionary);
        }
        final Iterator<VariantContext> posteriorsIter = reader.query(currentContig, 1, contigRecord.getSequenceLength());
        final IntervalTree<Map<String,double[]>> tree = new IntervalTree<>();
        while (posteriorsIter.hasNext()) {
            final VariantContext variant = posteriorsIter.next();
            final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd() - 1);
            final Map<String,double[]> samplePosteriorMap = new HashMap<>(SVUtils.hashMapCapacity(samples.size()));
            for (final Genotype genotype : variant.getGenotypes()) {
                final IntervalCopyNumberGenotypingData data = parseGenotype(genotype, interval);
                final double[] posteriors = new double[numCopyStates];
                int i = 0;
                for (final IntegerCopyNumberState state : copyStates) {
                    posteriors[i] = data.getCopyNumberPosteriorDistribution().getCopyNumberPosterior(state);
                    i++;
                }
                samplePosteriorMap.put(genotype.getSampleName(), posteriors);
            }
            tree.put(interval.getStart(), interval.getEnd(), samplePosteriorMap);
        }
        return tree;
    }

    private IntervalCopyNumberGenotypingData parseGenotype(final Genotype genotype, final SimpleInterval interval) {
        if (!genotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY)) {
            throw new UserException.BadInput("Copy number genotype not found in genotype: " + genotype.toString());
        }
        final String[] copyStatePLStrings = ((String)genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY)).split(",");

        //Posteriors reported as integer phred-scaled likelihoods and need to be renormalized
        final double[] approximatePosteriors = new double[copyStatePLStrings.length];
        double total = 0;
        for (int i = 0; i < copyStatePLStrings.length; i++) {
            final int likelihood = Integer.valueOf(copyStatePLStrings[i]);
            approximatePosteriors[i] = QualityUtils.qualToErrorProb(likelihood);
            total += approximatePosteriors[i];
        }

        final Map<IntegerCopyNumberState,Double> copyNumberPosteriors = new HashMap<>(SVUtils.hashMapCapacity(copyStatePLStrings.length));
        for (int i = 0; i < copyStatePLStrings.length; i++) {
            copyNumberPosteriors.put(new IntegerCopyNumberState(i), FastMath.log(Math.max(approximatePosteriors[i] / total, Double.MIN_VALUE)));
        }

        final String sample = genotype.getSampleName();
        return new IntervalCopyNumberGenotypingData(
                interval,
                new CopyNumberPosteriorDistribution(copyNumberPosteriors),
                getNeutralCopyState(sample, interval.getContig()));
    }

    private IntegerCopyNumberState getNeutralCopyState(final String sample, final String contig) {
        if (!sampleContigPloidyMap.containsKey(sample)) {
            throw new IllegalArgumentException("Sample " + sample + " not found in ploidy contig calls.");
        }
        final Map<String,Integer> contigPloidyMap = sampleContigPloidyMap.get(sample);
        if (!contigPloidyMap.containsKey(contig)) {
            throw new IllegalArgumentException("Contig '" + contig + "' not found in ploidy contig calls for sample '" + sample + "'.");
        }
        return new IntegerCopyNumberState(contigPloidyMap.get(contig));
    }

    private static Map<String,Integer> getContigToPloidyCallMap(final CalledContigPloidyCollection contigPloidyCollection) {
        return contigPloidyCollection.getRecords().stream().collect(Collectors.toMap(p -> p.getContig(), p -> p.getPloidy()));
    }

    public static Map<String,Map<String,Integer>> getSampleContigPloidyMap(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        return contigPloidyCollections.stream().collect(Collectors.toMap(p -> p.getMetadata().getSampleName(), p -> getContigToPloidyCallMap(p)));
    }
}
