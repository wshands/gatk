package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.*;
import java.util.stream.Collectors;

public abstract class SVCollapser<T extends SVCallRecord> {

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

    public static final BreakpointSummaryStrategy DEFAULT_STRATEGY = BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;
    private final BreakpointSummaryStrategy breakpointSummaryStrategy;

    public SVCollapser(final BreakpointSummaryStrategy breakpointSummaryStrategy) {
        this.breakpointSummaryStrategy = breakpointSummaryStrategy;
    }

    public SVCollapser() {
        this(DEFAULT_STRATEGY);
    }

    /**
     * Find a single base-class record representative of all the calls in a cluster of {@param items}.
     * @param items   the events that are clustered together
     * @return  a call approximating the average event for the cluster and containing all the algorithms and genotypes
     */
    abstract public T collapse(final Collection<T> items);
    abstract protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes);

    /**
     * Subclasses should call this from their implementations of {@link SVCollapser#collapse}.
     * @param items
     * @return
     */
    protected final SVCallRecord collapseToBaseRecord(final Collection<T> items) {
        final Collection<T> mostPreciseCalls = getMostPreciseCalls(items);
        final T exampleCall = mostPreciseCalls.iterator().next();
        final String id = collapseIds(items);
        final List<String> algorithms = collapseAlgorithms(items);
        final List<Genotype> clusterSamples = collapseAllGenotypes(items);
        final StructuralVariantType type = collapseTypes(items);
        final Map.Entry<Integer, Integer> coordinates = collapseInterval(items);
        final int start = coordinates.getKey();
        final int end = coordinates.getValue();
        final int length = collapseLength(mostPreciseCalls, start, end, type);
        return new SVCallRecord(id, exampleCall.getContigA(), start, exampleCall.getStrandA(),
                exampleCall.getContigB(), end, exampleCall.getStrandB(), type, length, algorithms, clusterSamples);
    }

    private List<Genotype> collapseAllGenotypes(final Collection<T> items) {
        return items.stream()
                .map(T::getGenotypes)
                .flatMap(GenotypesContext::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName))
                .values()
                .stream()
                .map(this::collapseSampleGenotypes)
                .collect(Collectors.toList());
    }

    protected final int collapseLength(final Collection<T> items, final int newStart, final int newEnd,
                                       final StructuralVariantType newType) {
        final T exampleCall = items.iterator().next();
        if (!exampleCall.isIntrachromosomal()) {
            return -1;
        } else if (newType.equals(StructuralVariantType.INS)) {
            return (int) Math.round(items.stream().mapToInt(T::getLength).average().getAsDouble());
        } else {
            return newEnd - newStart + 1;
        }
    }

    protected String collapseIds(final Collection<T> records) {
        return records.iterator().next().getId();
    }

    protected Collection<T> getMostPreciseCalls(final Collection<T> items) {
        if (items.stream().allMatch(SVClusterEngine::isDepthOnlyCall)) {
            return items;
        } else {
            return items.stream().filter(call -> !SVClusterEngine.isDepthOnlyCall(call)).collect(Collectors.toList());
        }
    }

    /**
     * @param items
     * @return (key, value) entry of (start, end)
     */
    protected Map.Entry<Integer,Integer> collapseInterval(final Collection<T> items) {
        final T exampleCall = items.iterator().next();
        if (breakpointSummaryStrategy == null) {
            return new AbstractMap.SimpleImmutableEntry<>(exampleCall.getPositionA(), exampleCall.getPositionB());
        }
        final List<Integer> startPositions = items.stream().map(SVCallRecord::getPositionA).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = items.stream().map(SVCallRecord::getPositionB).sorted().collect(Collectors.toList());
        //use the mid value of the sorted list so the start and end represent real breakpoint observations
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean;
        } else {
            switch (breakpointSummaryStrategy) {
                case MEDIAN_START_MEDIAN_END:
                    newStart = medianStart;
                    newEnd = medianEnd;
                    break;
                case MIN_START_MAX_END:
                    newStart = startPositions.get(0);
                    newEnd = endPositions.get(endPositions.size() - 1);
                    break;
                case MAX_START_MIN_END:
                    newStart = startPositions.get(startPositions.size() - 1);
                    newEnd = endPositions.get(0);
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
        return new AbstractMap.SimpleImmutableEntry<>(newStart, newEnd);
    }

    protected StructuralVariantType collapseTypes(final Collection<T> records) {
        final Set<StructuralVariantType> types = records.stream().map(T::getType).collect(Collectors.toSet());
        if (types.size() == 1) {
            return types.iterator().next();
        }
        if (types.stream().allMatch(SVClusterEngine::isCnvType)) {
            return StructuralVariantType.CNV;
        }
        throw new IllegalArgumentException("Incompatible SV types found in cluster");
    }

    protected List<String> collapseAlgorithms(final Collection<T> records) {
        return records.stream()
                .map(T::getAlgorithms)
                .flatMap(Collection::stream)
                .distinct()
                .collect(Collectors.toList());
    }
}
