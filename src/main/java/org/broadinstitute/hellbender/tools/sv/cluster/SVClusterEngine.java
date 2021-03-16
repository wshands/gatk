package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.function.BiPredicate;
import java.util.function.Function;

public class SVClusterEngine<T extends SVCallRecord> extends LocatableClusterEngine<T> {

    protected final boolean enableCNV;  // DEL/DUP multi-allelic site clustering
    protected ClusteringParameters depthOnlyParams;
    protected ClusteringParameters mixedParams;
    protected ClusteringParameters evidenceParams;

    protected static final ClusteringParameters DEFAULT_DEPTH_ONLY_PARAMS =
            new DepthClusteringParameters(0.8, 0);
    protected static final ClusteringParameters DEFAULT_MIXED_PARAMS =
            new MixedClusteringParameters(0.8, 1000);
    protected static final ClusteringParameters DEFAULT_EVIDENCE_PARAMS =
            new EvidenceClusteringParameters(0.5, 500);

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
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        if (a.getType() != b.getType()) {
            if (!enableCNV) {
                // CNV clustering disabled, so no type mixing
                return false;
            } else if (!(isCnvType(a.getType()) && isCnvType(b.getType()))) {
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
        return clusterTogetherWithParams(a, b, evidenceParams)
                || clusterTogetherWithParams(a, b, depthOnlyParams)
                || clusterTogetherWithParams(a, b, mixedParams);
    }

    private boolean clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b, final ClusteringParameters params) {
        // Type check
        if (!params.isValidPair(a, b)) return false;

        // Reciprocal overlap
        final boolean isOverlap;
        if (a.isIntrachromosomal()) {
            final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionB());
            final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionB());
            isOverlap = IntervalUtils.isReciprocalOverlap(intervalA, intervalB, params.getReciprocalOverlap());
        } else {
            isOverlap = true;
        }

        // Breakend proximity
        final SimpleInterval intervalA1 = a.getPositionAInterval().expandWithinContig(params.getWindow(), dictionary);
        final SimpleInterval intervalB1 = b.getPositionAInterval().expandWithinContig(params.getWindow(), dictionary);
        final SimpleInterval intervalA2 = a.getPositionBInterval().expandWithinContig(params.getWindow(), dictionary);
        final SimpleInterval intervalB2 = b.getPositionBInterval().expandWithinContig(params.getWindow(), dictionary);
        final boolean isProximity = intervalA1.overlaps(intervalB1) && intervalA2.overlaps(intervalB2);
        if (params.isOverlapAndProximity()) {
            return isOverlap && isProximity;
        } else {
            return isOverlap || isProximity;
        }
    }

    @Override
    protected int getMaxClusterableStartingPosition(final SVCallRecord call) {
        final String contig = call.getContigA();
        final boolean isDepthOnly = call.isDepthOnlyCall();
        final int contigLength = dictionary.getSequence(contig).getSequenceLength();
        // Reciprocal overlap window
        final int maxPositionByOverlap;
        if (call.isIntrachromosomal()) {
            final double overlap = Math.min(isDepthOnly ? depthOnlyParams.getReciprocalOverlap() : evidenceParams.getReciprocalOverlap(), mixedParams.getReciprocalOverlap());
            final int maxPosition = (int) (call.getPositionA() + (1.0 - overlap) * (call.getPositionB() - call.getPositionA()));
            maxPositionByOverlap = Math.min(maxPosition, contigLength);
        } else {
            maxPositionByOverlap = call.getPositionA();
        }

        // Breakend proximity window
        final int window = Math.max(isDepthOnly ? depthOnlyParams.getWindow() : evidenceParams.getWindow(), mixedParams.getWindow());
        final int maxPositionByWindow = Math.min(call.getPositionA() + window, contigLength);

        if (isDepthOnly) {
            return Math.max(maxPositionByOverlap, maxPositionByWindow);
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

    public static boolean isCnvType(final StructuralVariantType type) {
        return type == StructuralVariantType.DEL || type == StructuralVariantType.DUP || type == StructuralVariantType.CNV;
    }

    public static class ClusteringParameters {

        private final double reciprocalOverlap;
        private final int window;
        private final boolean overlapAndProximity;
        private final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate;

        public ClusteringParameters(final double reciprocalOverlap, final int window, final boolean overlapAndProximity,
                                    final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate) {
            this.reciprocalOverlap = reciprocalOverlap;
            this.window = window;
            this.overlapAndProximity = overlapAndProximity;
            this.validRecordsPredicate = validRecordsPredicate;
        }

        public double getReciprocalOverlap() {
            return reciprocalOverlap;
        }

        public int getWindow() {
            return window;
        }

        public boolean isOverlapAndProximity() {
            return overlapAndProximity;
        }

        public boolean isValidPair(final SVCallRecord a, final SVCallRecord b) {
            return validRecordsPredicate.test(a, b);
        }
    }

    public static final class DepthClusteringParameters extends ClusteringParameters {
        public DepthClusteringParameters(final double reciprocalOverlap, final int window) {
            super(reciprocalOverlap, window, false, (a,b) -> a.isDepthOnlyCall() && b.isDepthOnlyCall());
        }
    }

    public static final class EvidenceClusteringParameters extends ClusteringParameters {
        public EvidenceClusteringParameters(final double reciprocalOverlap, final int window) {
            super(reciprocalOverlap, window, true, (a,b) -> !a.isDepthOnlyCall() && !b.isDepthOnlyCall());
        }
    }

    public static final class MixedClusteringParameters extends ClusteringParameters {
        public MixedClusteringParameters(final double reciprocalOverlap, final int window) {
            super(reciprocalOverlap, window, true, (a,b) -> a.isDepthOnlyCall() != b.isDepthOnlyCall());
        }
    }

}
