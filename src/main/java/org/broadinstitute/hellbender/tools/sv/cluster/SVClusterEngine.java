package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.COPY_NUMBER_FORMAT;

public class SVClusterEngine<T extends SVCallRecord> extends LocatableClusterEngine<T> {

    protected final boolean enableCNV;  // DEL/DUP multi-allelic site clustering
    protected ClusteringParameters depthOnlyParams;
    protected ClusteringParameters mixedParams;
    protected ClusteringParameters evidenceParams;

    public static final int INSERTION_ASSUMED_LENGTH_FOR_OVERLAP = 50;

    public static final ClusteringParameters DEFAULT_DEPTH_ONLY_PARAMS =
            new DepthClusteringParameters(0.8, 0, 0);
    public static final ClusteringParameters DEFAULT_MIXED_PARAMS =
            new MixedClusteringParameters(0.8, 1000, 0);
    public static final ClusteringParameters DEFAULT_EVIDENCE_PARAMS =
            new EvidenceClusteringParameters(0.5, 500, 0);

    public SVClusterEngine(final SAMSequenceDictionary dictionary,
                           final CLUSTERING_TYPE clusteringType,
                           final boolean enableCNV,
                           final Function<Collection<T>, T> collapser) {
        super(dictionary, clusteringType, collapser);
        this.depthOnlyParams = DEFAULT_DEPTH_ONLY_PARAMS;
        this.mixedParams = DEFAULT_MIXED_PARAMS;
        this.evidenceParams = DEFAULT_EVIDENCE_PARAMS;
        this.enableCNV = enableCNV;
    }

    @Override
    boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        if (a.getType() != b.getType()) {
            if (!enableCNV) {
                // CNV clustering disabled, so no type mixing
                return false;
            } else if (!(a.isCNV() && b.isCNV())) {
                // CNV clustering enabled, but at least one was not a CNV type
                return false;
            }
        }
        // Contigs match
        if (!(a.getContigA().equals(b.getContigA()) && a.getContigB().equals(b.getContigB()))) {
            return false;
        }
        // Strands match
        if (a.getStrandA() != b.getStrandA() || a.getStrandB() != b.getStrandB()) {
            return false;
        }
        // Checks appropriate parameter set
        return clusterTogetherWithParams(a, b, evidenceParams, dictionary)
                || clusterTogetherWithParams(a, b, depthOnlyParams, dictionary)
                || clusterTogetherWithParams(a, b, mixedParams, dictionary);
    }

    private static Collection<Genotype> getCNVCarrierGenotypesByCopyNumber(final SVCallRecord record) {
        Utils.validate(record.isCNV(), "Cannot determine carriers of non-CNV using copy number attribute.");
        final List<Allele> altAlleles = record.getAltAlleles();
        Utils.validate(altAlleles.size() <= 1,
                "Carrier samples cannot be determined by copy number for multi-allelic sites. Set sample overlap threshold to 0.");
        if (altAlleles.isEmpty()) {
            return Collections.emptyList();
        }
        final Allele altAllele = altAlleles.get(0);
        return record.getGenotypes().stream()
                .filter(g -> isCNVCarrierByCopyNumber(g, altAllele))
                .collect(Collectors.toList());
    }

    private static boolean isCNVCarrierByCopyNumber(final Genotype genotype, final Allele altAllele) {
        final int ploidy = genotype.getPloidy();
        if (ploidy == 0 || !genotype.hasExtendedAttribute(COPY_NUMBER_FORMAT)) {
            return false;
        }
        final int copyNumber = VariantContextGetters.getAttributeAsInt(genotype, COPY_NUMBER_FORMAT, 0);
        if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
            return copyNumber < ploidy;
        } else {
            // DUP
            return copyNumber > ploidy;
        }
    }

    private static Set<String> getCarrierSamplesByCopyNumber(final SVCallRecord record) {
        return getCNVCarrierGenotypesByCopyNumber(record).stream().map(Genotype::getSampleName).collect(Collectors.toSet());
    }

    private static boolean hasSampleSetOverlap(final Set<String> samplesA, final Set<String> samplesB, final double minSampleOverlap) {
        final int denom = Math.max(samplesA.size(), samplesB.size());
        if (denom == 0) {
            return true;
        }
        final double sampleOverlap = getSampleSetOverlap(samplesA, samplesB) / (double) denom;
        return sampleOverlap >= minSampleOverlap;
    }

    private static int getSampleSetOverlap(final Collection<String> a, final Set<String> b) {
        return (int) a.stream().filter(b::contains).count();
    }

    private static boolean hasDefinedCopyNumbers(final Collection<Genotype> genotypes) {
        return genotypes.stream().anyMatch(g -> g.hasExtendedAttribute(COPY_NUMBER_FORMAT));
    }

    private static boolean hasExplicitGenotypes(final Collection<Genotype> genotypes) {
        return genotypes.stream().anyMatch(Genotype::isCalled);
    }

    private static Set<String> getCarrierSamplesByGenotype(final SVCallRecord record) {
        final Set<Allele> altAlleles = new HashSet<>(record.getAltAlleles());
        return record.getGenotypes().stream()
                .filter(g -> g.getAlleles().stream().anyMatch(altAlleles::contains))
                .map(Genotype::getSampleName)
                .collect(Collectors.toSet());
    }

    protected static Set<String> getCarrierSamples(final SVCallRecord record) {
        if (record.isCNV() && !hasExplicitGenotypes(record.getGenotypes()) && hasDefinedCopyNumbers(record.getGenotypes())) {
            return getCarrierSamplesByCopyNumber(record);
        } else {
            return getCarrierSamplesByGenotype(record);
        }
    }

    protected static boolean hasSampleOverlap(final SVCallRecord a, final SVCallRecord b, final double minSampleOverlap) {
        if (minSampleOverlap > 0) {
            final Set<String> samplesA = getCarrierSamples(a);
            final Set<String> samplesB = getCarrierSamples(b);
            return hasSampleSetOverlap(samplesA, samplesB, minSampleOverlap);
        } else {
            return true;
        }
    }

    private static boolean clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b,
                                                     final ClusteringParameters params,
                                                     final SAMSequenceDictionary dictionary) {
        // Type check
        if (!params.isValidPair(a, b)) {
            return false;
        }

        // Sample overlap
        if (!hasSampleOverlap(a, b, params.getSampleOverlap())) {
            return false;
        }

        // Reciprocal overlap
        final boolean isOverlap;
        if (a.isIntrachromosomal()) {
            final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionA() + getLengthForOverlap(a) - 1);
            final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionA() + getLengthForOverlap(b) - 1);
            isOverlap = IntervalUtils.isReciprocalOverlap(intervalA, intervalB, params.getReciprocalOverlap());
        } else {
            isOverlap = true;
        }
        if (params.isOverlapAndProximity() && !isOverlap) {
            return false;
        } else if (!params.isOverlapAndProximity() && isOverlap) {
            return true;
        }

        // Breakend proximity
        final SimpleInterval intervalA1 = a.getPositionAInterval().expandWithinContig(params.getWindow(), dictionary);
        final SimpleInterval intervalA2 = a.getPositionBInterval().expandWithinContig(params.getWindow(), dictionary);
        final SimpleInterval intervalB1 = b.getPositionAInterval();
        final SimpleInterval intervalB2 = b.getPositionBInterval();
        return intervalA1.overlaps(intervalB1) && intervalA2.overlaps(intervalB2);
    }

    private static int getLengthForOverlap(final SVCallRecord record) {
        Utils.validate(record.isIntrachromosomal(), "Record even must be intra-chromosomal");
        if (record.getType() == StructuralVariantType.INS) {
            return record.getLength() < 1 ? INSERTION_ASSUMED_LENGTH_FOR_OVERLAP : record.getLength();
        } else {
            // TODO lengths less than 1 shouldn't be valid
            return Math.max(record.getLength(), 1);
        }
    }

    @Override
    protected int getMaxClusterableStartingPosition(final SVCallRecord call) {
        return Math.max(
                getMaxClusterableStartingPositionWithParams(call, call.isDepthOnly() ? depthOnlyParams : evidenceParams, dictionary),
                getMaxClusterableStartingPositionWithParams(call, mixedParams, dictionary)
        );
    }

    private static int getMaxClusterableStartingPositionWithParams(final SVCallRecord call, final ClusteringParameters params,
                                                            final SAMSequenceDictionary dictionary) {
        final String contig = call.getContigA();
        final int contigLength = dictionary.getSequence(contig).getSequenceLength();
        // Reciprocal overlap window
        final int maxPositionByOverlap;
        if (call.isIntrachromosomal()) {
            final int maxPosition = (int) (call.getPositionA() + (1.0 - params.getReciprocalOverlap()) * getLengthForOverlap(call));
            maxPositionByOverlap = Math.min(maxPosition, contigLength);
        } else {
            maxPositionByOverlap = call.getPositionA();
        }

        // Breakend proximity window
        final int maxPositionByWindow = Math.min(call.getPositionA() + params.getWindow(), contigLength);

        if (params.isOverlapAndProximity()) {
            return Math.min(maxPositionByOverlap, maxPositionByWindow);
        } else {
            return Math.max(maxPositionByOverlap, maxPositionByWindow);
        }
    }

    public final ClusteringParameters getDepthOnlyParams() {
        return depthOnlyParams;
    }
    public final ClusteringParameters getMixedParams() {
        return mixedParams;
    }
    public final ClusteringParameters getEvidenceParams() {
        return evidenceParams;
    }

    public final void setDepthOnlyParams(ClusteringParameters depthOnlyParams) { this.depthOnlyParams = depthOnlyParams; }
    public final void setMixedParams(ClusteringParameters mixedParams) {
        this.mixedParams = mixedParams;
    }
    public final void setEvidenceParams(ClusteringParameters evidenceParams) {
        this.evidenceParams = evidenceParams;
    }

    public static class ClusteringParameters {

        private final double reciprocalOverlap;
        private final int window;
        private final double sampleOverlap;
        private final boolean overlapAndProximity;
        private final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate;

        public ClusteringParameters(final double reciprocalOverlap, final int window, final double sampleOverlap,
                                    final boolean overlapAndProximity, final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate) {
            this.reciprocalOverlap = reciprocalOverlap;
            this.window = window;
            this.sampleOverlap = sampleOverlap;
            this.overlapAndProximity = overlapAndProximity;
            this.validRecordsPredicate = validRecordsPredicate;
        }

        public double getReciprocalOverlap() {
            return reciprocalOverlap;
        }

        public int getWindow() {
            return window;
        }

        public double getSampleOverlap() {
            return sampleOverlap;
        }

        public boolean isOverlapAndProximity() {
            return overlapAndProximity;
        }

        public boolean isValidPair(final SVCallRecord a, final SVCallRecord b) {
            return validRecordsPredicate.test(a, b);
        }
    }

    public static final class DepthClusteringParameters extends ClusteringParameters {
        public DepthClusteringParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
            super(reciprocalOverlap, window, sampleOverlap, false, (a,b) -> a.isDepthOnly() && b.isDepthOnly());
        }
    }

    public static final class EvidenceClusteringParameters extends ClusteringParameters {
        public EvidenceClusteringParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
            super(reciprocalOverlap, window, sampleOverlap, true, (a,b) -> !a.isDepthOnly() && !b.isDepthOnly());
        }
    }

    public static final class MixedClusteringParameters extends ClusteringParameters {
        public MixedClusteringParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
            super(reciprocalOverlap, window, sampleOverlap, true, (a,b) -> a.isDepthOnly() != b.isDepthOnly());
        }
    }

}
