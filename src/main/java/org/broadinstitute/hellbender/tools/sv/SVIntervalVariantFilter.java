package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Keep only variants meeting specified overlap requirements with the input intervals.
 */
public final class SVIntervalVariantFilter implements VariantFilter {

    private final Map<String, IntervalTree<Object>> includedIntervalsTreeMap;
    private final boolean requireBreakendOverlap;
    private final double minOverlapFraction;

    // Command line parser requires a no-arg constructor
    public SVIntervalVariantFilter() {
        this.includedIntervalsTreeMap = null;
        this.requireBreakendOverlap = false;
        this.minOverlapFraction = 0;
    }

    public SVIntervalVariantFilter(final Collection<SimpleInterval> intervals,
                                   final boolean requireBreakendOverlap,
                                   final double minOverlapFraction) {
        Utils.validateArg(0 <= minOverlapFraction && minOverlapFraction <= 1, "Invalid overlap fraction: " + minOverlapFraction);
        this.requireBreakendOverlap = requireBreakendOverlap;
        this.minOverlapFraction = minOverlapFraction;
        includedIntervalsTreeMap = new HashMap<>();
        for (final SimpleInterval interval : intervals) {
            includedIntervalsTreeMap.putIfAbsent(interval.getContig(), new IntervalTree<>());
            includedIntervalsTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
        }
    }

    @Override
    public boolean test(final VariantContext variant) {
        final String contig1 = variant.getContig();
        final int start = variant.getStart();
        final IntervalTree<Object> startTree = includedIntervalsTreeMap.get(contig1);
        // First contig included
        if (startTree == null) {
            return false;
        }
        final String contig2 = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, contig1);
        final IntervalTree<Object> endTree = includedIntervalsTreeMap.get(contig2);
        // Second contig included
        if (endTree == null) {
            return false;
        }
        // Breakends both included, if required
        final int end = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, variant.getEnd());
        if (requireBreakendOverlap && !(startTree.overlappers(start, start + 1).hasNext()
                && endTree.overlappers(end, end + 1).hasNext())) {
            return false;
        }
        // Can include all interchromosomal variants at this point
        if (!variant.getContig().equals(contig2)) {
            return true;
        }
        // Check overlap fraction
        final long overlapLength = totalOverlap(start, end, startTree);
        final double overlapFraction = overlapLength / (double) (end - variant.getStart());
        return overlapFraction >= minOverlapFraction;
    }

    protected static <T> long totalOverlap(final int start, final int end, final IntervalTree<T> tree) {
        final Iterator<IntervalTree.Node<T>> iter = tree.overlappers(start, end);
        long overlap = 0;
        while (iter.hasNext()) {
            final IntervalTree.Node<T> node = iter.next();
            overlap += intersectionLength(start, end, node.getStart(), node.getEnd());
        }
        return overlap;
    }

    protected static long intersectionLength(final int start1, final int end1, final int start2, final int end2) {
        return Math.max(0, Math.min(end1, end2) - Math.max(start1, start2) + 1);
    }
}