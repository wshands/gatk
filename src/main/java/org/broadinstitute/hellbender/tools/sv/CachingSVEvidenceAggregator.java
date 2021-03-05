package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.utils.*;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public abstract class CachingSVEvidenceAggregator<T extends Feature> {

    private final FeatureDataSource<T> source;
    private SimpleInterval cacheInterval;
    protected final SAMSequenceDictionary dictionary;
    protected final ProgressMeter progressMeter;

    public CachingSVEvidenceAggregator(final FeatureDataSource<T> source,
                                       final SAMSequenceDictionary dictionary,
                                       final ProgressMeter progressMeter) {
        this.source = source;
        this.dictionary = dictionary;
        this.cacheInterval = null;
        this.progressMeter = progressMeter;
    }

    abstract protected SimpleInterval getEvidenceQueryInterval(final SVCallRecordWithEvidence record);
    abstract protected SVCallRecordWithEvidence assignEvidence(final SVCallRecordWithEvidence call, final List<T> evidence);
    protected boolean evidenceForCallFilter(final SVCallRecord record, final T evidence) { return true; }
    protected boolean evidenceFilter(final T evidence) { return true; }

    public Stream<SVCallRecordWithEvidence> collectEvidence(final List<SVCallRecordWithEvidence> calls) {
        Utils.nonNull(calls);
        final OverlapDetector<SimpleInterval> overlapDetector = getEvidenceOverlapDetector(calls);
        return calls.stream()
                .collect(Collectors.toList())
                .stream()
                .map(r -> collectRecordEvidence(r, overlapDetector));
    }

    private final SVCallRecordWithEvidence collectRecordEvidence(final SVCallRecordWithEvidence record,
                                                                 final OverlapDetector<SimpleInterval> overlapDetector) {
        final List<T> evidence = getEvidenceOnInterval(getEvidenceQueryInterval(record), overlapDetector).stream()
                .filter(e -> evidenceForCallFilter(record, e))
                .collect(Collectors.toList());
        return assignEvidence(record, evidence);
    }

    private final List<T> getEvidenceOnInterval(final SimpleInterval interval,
                                          final OverlapDetector<SimpleInterval> overlapDetector) {
        if (invalidCacheInterval(cacheInterval, interval)) {
            Utils.nonNull(overlapDetector, "Evidence cache missed but overlap detector is null");
            final Set<SimpleInterval> queryIntervalSet = overlapDetector.getOverlaps(interval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Dvidence interval " + interval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            cacheInterval = queryIntervalSet.iterator().next();
            source.queryAndPrefetch(cacheInterval);
        }
        return source.queryAndPrefetch(interval).stream()
                .filter(this::evidenceFilter)
                .collect(Collectors.toList());
    }

    private final OverlapDetector<SimpleInterval> getEvidenceOverlapDetector(final List<SVCallRecordWithEvidence> calls) {
        final List<SimpleInterval> rawIntervals = calls.stream()
                .map(this::getEvidenceQueryInterval)
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .collect(Collectors.toList());
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> rawLocs = IntervalUtils.genomeLocsFromLocatables(parser, rawIntervals);
        final List<GenomeLoc> mergedLocs = IntervalUtils.mergeIntervalLocations(rawLocs, IntervalMergingRule.ALL);
        final List<SimpleInterval> mergedIntervals = IntervalUtils.convertGenomeLocsToSimpleIntervals(mergedLocs);
        return OverlapDetector.create(mergedIntervals);
    }

    private final boolean invalidCacheInterval(final SimpleInterval cacheInterval, final SimpleInterval queryInterval) {
        return cacheInterval == null
                || !queryInterval.getContig().equals(cacheInterval.getContig())
                || !queryInterval.spanWith(cacheInterval).equals(cacheInterval);
    }
}
