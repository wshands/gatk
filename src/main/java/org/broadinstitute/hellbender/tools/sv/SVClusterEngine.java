package org.broadinstitute.hellbender.tools.sv;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SVClusterEngine extends LocatableClusterEngine<SVCallRecord> {

    public enum BreakpointSummaryStrategy {
        /**
         * Use the (first) middle value to summarize cluster starts and ends, such that the start and end were seen in the data
         */
        MEDIAN_START_MEDIAN_END,

        /**
         * A conservative strategy to summarize a cluster by its smallest extent
         */
        MIN_START_MAX_END,

        /**
         * A permissive strategy to summarize a cluster by it largest extent
         */
        MAX_START_MIN_END,

        /**
         * Summarize a cluster using the mean value for each end, even if that value was not represented in any sample
         */
        MEAN_START_MEAN_END

    }

    protected static class ClusteringParameters {

        private final double reciprocalOverlap;
        private final int window;
        private final int padding;
        private final boolean overlapAndProximity;
        private final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate;

        public ClusteringParameters(final double reciprocalOverlap, final int window, final int padding,
                                    final boolean overlapAndProximity,
                                    final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate) {
            this.reciprocalOverlap = reciprocalOverlap;
            this.window = window;
            this.padding = padding;
            this.overlapAndProximity = overlapAndProximity;
            this.validRecordsPredicate = validRecordsPredicate;
        }

        public double getReciprocalOverlap() {
            return reciprocalOverlap;
        }

        public int getWindow() {
            return window;
        }

        public int getPadding() {
            return padding;
        }

        public boolean isOverlapAndProximity() {
            return overlapAndProximity;
        }

        public boolean isValidPair(final SVCallRecord a, final SVCallRecord b) {
            return validRecordsPredicate.test(a, b);
        }
    }

    public static class DepthClusteringParameters extends ClusteringParameters {
        public DepthClusteringParameters(final double reciprocalOverlap, final int window, final int padding) {
            super(reciprocalOverlap, window, padding, false, (a,b) -> isDepthOnlyCall(a) && isDepthOnlyCall(b));
        }
    }

    public static class EvidenceClusteringParameters extends ClusteringParameters {
        public EvidenceClusteringParameters(final double reciprocalOverlap, final int window, final int padding) {
            super(reciprocalOverlap, window, padding, true, (a,b) -> !isDepthOnlyCall(a) && !isDepthOnlyCall(b));
        }
    }

    public static class MixedClusteringParameters extends ClusteringParameters {
        public MixedClusteringParameters(final double reciprocalOverlap, final int window, final int padding) {
            super(reciprocalOverlap, window, padding, true, (a,b) -> isDepthOnlyCall(a) != isDepthOnlyCall(b));
        }
    }

    protected static final ClusteringParameters DEFAULT_DEPTH_ONLY_PARAMS =
            new DepthClusteringParameters(0.8, 1000, 0);
    protected static final ClusteringParameters DEFAULT_MIXED_PARAMS =
            new MixedClusteringParameters(0.5, 1000, 1000);
    protected static final ClusteringParameters DEFAULT_EVIDENCE_PARAMS =
            new EvidenceClusteringParameters(0.5, 500, 500);

    protected ClusteringParameters depthOnlyParams;
    protected ClusteringParameters mixedParams;
    protected ClusteringParameters evidenceParams;
    protected BreakpointSummaryStrategy breakpointSummaryStrategy;

    public SVClusterEngine(final SAMSequenceDictionary dictionary) {
        super(dictionary, CLUSTERING_TYPE.MAX_CLIQUE, null);
        this.breakpointSummaryStrategy = BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;
        this.depthOnlyParams = DEFAULT_DEPTH_ONLY_PARAMS;
        this.mixedParams = DEFAULT_MIXED_PARAMS;
        this.evidenceParams = DEFAULT_EVIDENCE_PARAMS;
    }

    public SVClusterEngine(final SAMSequenceDictionary dictionary, boolean singleLinkage, BreakpointSummaryStrategy strategy) {
        super(dictionary, singleLinkage ? CLUSTERING_TYPE.SINGLE_LINKAGE : CLUSTERING_TYPE.MAX_CLIQUE, null);
        this.breakpointSummaryStrategy = strategy;
    }

    public BreakpointSummaryStrategy getBreakpointSummaryStrategy() {
        return breakpointSummaryStrategy;
    }

    public void setBreakpointSummaryStrategy(BreakpointSummaryStrategy breakpointSummaryStrategy) {
        this.breakpointSummaryStrategy = breakpointSummaryStrategy;
    }

    public ClusteringParameters getDepthOnlyParams() {
        return depthOnlyParams;
    }

    public void setDepthOnlyParams(ClusteringParameters depthOnlyParams) {
        this.depthOnlyParams = depthOnlyParams;
    }

    public ClusteringParameters getMixedParams() {
        return mixedParams;
    }

    public void setMixedParams(ClusteringParameters mixedParams) {
        this.mixedParams = mixedParams;
    }

    public ClusteringParameters getEvidenceParams() {
        return evidenceParams;
    }

    public void setEvidenceParams(ClusteringParameters evidenceParams) {
        this.evidenceParams = evidenceParams;
    }

    /**
     * Find a single call representative of all the calls in the {@param cluster}
     * @param cluster   the events that are clustered together
     * @return  a call approximating the average event for the cluster and containing all the algorithms and genotypes
     */

    @Override
    protected SVCallRecord flattenCluster(final Collection<SVCallRecord> cluster) {
        final Collection<SVCallRecord> mostPreciseCalls;
        if (cluster.stream().allMatch(SVClusterEngine::isDepthOnlyCall)) {
            mostPreciseCalls = cluster;
        } else {
            mostPreciseCalls = cluster.stream().filter(call -> !SVDepthOnlyCallDefragmenter.isDepthOnlyCall(call)).collect(Collectors.toList());
        }
        final List<Integer> startPositions = mostPreciseCalls.stream().map(SVCallRecord::getPositionA).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = mostPreciseCalls.stream().map(SVCallRecord::getPositionB).sorted().collect(Collectors.toList());
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final SVCallRecord exampleCall = mostPreciseCalls.iterator().next();
        final int length = exampleCall.getContigA().equals(exampleCall.getContigB()) && !exampleCall.getType().equals(StructuralVariantType.INS) ? medianEnd - medianStart + 1 : exampleCall.getLength();
        final List<String> algorithms = cluster.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final List<Genotype> clusterSamples = cluster.stream().flatMap(v -> v.getGenotypes().stream()).collect(Collectors.toList());
        final List<StructuralVariantType> observedTypes = cluster.stream().map(SVCallRecord::getType).distinct().collect(Collectors.toList());
        final StructuralVariantType clusterType = observedTypes.size() == 1 ? observedTypes.get(0) : StructuralVariantType.CNV;

        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean;
        } else {
            switch (this.breakpointSummaryStrategy) {
                case MEDIAN_START_MEDIAN_END:
                    newStart = medianStart;
                    newEnd = medianEnd;
                    break;
                case MIN_START_MAX_END:
                    newStart = startPositions.stream().min(Integer::compareTo).orElse(startPositions.get(0));
                    newEnd = endPositions.stream().max(Integer::compareTo).orElse(endPositions.get(0));
                    break;
                case MAX_START_MIN_END:
                    newStart = startPositions.stream().max(Integer::compareTo).orElse(startPositions.get(0));
                    newEnd = endPositions.stream().min(Integer::compareTo).orElse(endPositions.get(0));
                    break;
                case MEAN_START_MEAN_END:
                    newStart = (int)Math.round(new Mean().evaluate(Doubles.toArray(startPositions)));
                    newEnd = (int)Math.round(new Mean().evaluate(Doubles.toArray(endPositions)));
                    break;
                default:
                    newStart = medianStart;
                    newEnd = medianEnd;
            }

        }

        //TODO: merge evidence for WGS data
        return new SVCallRecord(exampleCall.getId(), exampleCall.getContigA(), newStart, exampleCall.getStrandA(),
                exampleCall.getContigB(), newEnd, exampleCall.getStrandB(), clusterType, length, algorithms, clusterSamples);
    }

    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        //if (!a.getType().equals(b.getType())) return false;  //TODO: do we need to keep dels and dupes separate?
        return clusterTogetherWithParams(a, b, evidenceParams)
                || clusterTogetherWithParams(a, b, depthOnlyParams)
                || clusterTogetherWithParams(a, b, mixedParams);
    }

    protected boolean clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b, final ClusteringParameters params) {
        // Type check
        if (!params.isValidPair(a, b)) return false;

        // Contigs match
        if (!(a.getContigA().equals(b.getContigA()) && a.getContigB().equals(b.getContigB()))) return false;

        // Reciprocal overlap
        final boolean isOverlap;
        if (a.isIntrachromosomal()) {
            final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionB()).expandWithinContig(params.getPadding(), dictionary);
            final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionB()).expandWithinContig(params.getPadding(), dictionary);
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

    /**
     * Determine an overlap interval for clustering using reciprocal overlap or breakend window, as applicable
     * Returned interval represents the interval in which the start position of a new event must fall in order to be added to the cluster (including {@param call})
     * @param call  new event to be clustered
     * @param clusterMinStartInterval    the cluster of interest, may be null
     * @return  an interval describing the cluster after {@param call} is added
     */
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecord call, final SimpleInterval clusterMinStartInterval) {
        final int minStart;
        final int maxStart;
        if (isDepthOnlyCall(call)) {
            minStart = (int) (call.getPositionB() - call.getLength() / depthOnlyParams.getReciprocalOverlap()); //start of an overlapping event such that call represents (reciprocal overlap) of that event
            maxStart = (int) (call.getPositionA() + (1.0 - depthOnlyParams.getReciprocalOverlap()) * call.getLength());
        } else {
            final int window = Math.max(mixedParams.getPadding(), Math.max(mixedParams.getWindow(), Math.max(evidenceParams.getWindow(), evidenceParams.getPadding())));
            minStart = call.getPositionA() - window;
            maxStart = call.getPositionA() + window;
        }
        final String currentContig = getCurrentContig();
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        //NOTE: this is an approximation -- best method would back calculate cluster bounds, then rederive start and end based on call + cluster
        final int newMinStart = Math.min(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.max(maxStart, clusterMinStartInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxStart, dictionary.getSequence(currentContig).getSequenceLength());
    }

    @Override
    protected SVDeduplicator<SVCallRecord> getDeduplicator() {
        final Function<Collection<SVCallRecord>, SVCallRecord> collapser = items -> SVCallRecordUtils.deduplicateWithRawCallAttribute(items, SVCallRecordUtils.ALLELE_COLLAPSER_DIPLOID_NO_CALL);
        return new SVCallRecordDeduplicator<>(collapser, dictionary);
    }

    public static boolean isDepthOnlyCall(final SVCallRecord call) {
        return call.getAlgorithms().size() == 1 && call.getAlgorithms().get(0).equals(GATKSVVCFConstants.DEPTH_ALGORITHM);
    }

    public double getMinReciprocalOverlap() {
        return depthOnlyParams.getReciprocalOverlap();
    }
}
