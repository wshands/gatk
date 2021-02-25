package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class SVClusterEngineNoCNV extends SVClusterEngine {

    public SVClusterEngineNoCNV(final SAMSequenceDictionary dictionary) {
        super(dictionary);
    }

    public SVClusterEngineNoCNV(final SAMSequenceDictionary dictionary, boolean singleLinkage, BreakpointSummaryStrategy strategy) {
        super(dictionary, singleLinkage, strategy);
    }

    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        if (!a.getType().equals(b.getType()) || a.getStrandA() != b.getStrandA() || a.getStrandB() != b.getStrandB()) {
            return false;
        }
        return super.clusterTogether(a, b);
    }

    // TODO optimize intervals
    @Override
    protected SimpleInterval getClusteringInterval(final SVCallRecord call, final SimpleInterval clusterMinStartInterval) {
        final SimpleInterval itemInterval = getClusteringInterval(call);
        if (clusterMinStartInterval == null) {
            return itemInterval;
        } else {
            return restrictIntervalsForClustering(itemInterval, clusterMinStartInterval);
        }
    }

    protected SimpleInterval restrictIntervalsForClustering(final SimpleInterval a, final SimpleInterval b) {
        Utils.validateArg(a.getContig().equals(b.getContig()), "Intervals are on different contigs");
        final int start;
        final int end;
        if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
            start = Math.max(a.getStart(), b.getStart());
            end = Math.min(a.getEnd(), b.getEnd());
        } else {
            start = Math.min(a.getStart(), b.getStart());
            end = Math.max(a.getEnd(), b.getEnd());
        }
        final String contig = a.getContig();
        return IntervalUtils.trimIntervalToContig(contig, start, end, dictionary.getSequence(contig).getSequenceLength());
    }

    protected SimpleInterval getClusteringInterval(final SVCallRecord call) {
        final String contig = call.getContigA();
        final boolean isDepthOnly = isDepthOnlyCall(call);
        // Reciprocal overlap window
        final SimpleInterval overlapInterval;
        if (call.isIntrachromosomal()) {
            final int padding = Math.max(isDepthOnly ? depthOnlyParams.getPadding() : evidenceParams.getPadding(), mixedParams.getPadding());
            final double overlap = Math.min(isDepthOnly ? depthOnlyParams.getReciprocalOverlap() : evidenceParams.getReciprocalOverlap(), mixedParams.getReciprocalOverlap());
            final SimpleInterval spanningInterval = new SimpleInterval(contig, call.getPositionA(), call.getPositionB()).expandWithinContig(padding, dictionary);
            final int start = (int) (spanningInterval.getEnd() - (spanningInterval.getLengthOnReference() / overlap));
            final int end = (int) (spanningInterval.getStart() + (1.0 - overlap) * spanningInterval.getLengthOnReference());
            overlapInterval = IntervalUtils.trimIntervalToContig(contig, start, end, dictionary.getSequence(contig).getSequenceLength());
        } else {
            overlapInterval = call.getPositionAInterval();
        }

        // Breakend proximity window
        final int window = Math.max(isDepthOnly ? depthOnlyParams.getWindow() : evidenceParams.getWindow(), mixedParams.getWindow());
        final SimpleInterval breakendInterval = call.getPositionAInterval().expandWithinContig(window, dictionary);

        return restrictIntervalsForClustering(overlapInterval, breakendInterval);
    }
}
