package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.*;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class BinnedCNVDefragmenter extends CNVDefragmenter {

    protected final TreeMap<GenomeLoc, Integer> genomicToBinMap;
    protected final List<GenomeLoc> coverageIntervals;
    final GenomeLocParser parser;

    //for single-sample clustering case
    public BinnedCNVDefragmenter(final SAMSequenceDictionary dictionary, final double paddingFraction,
                                 final double minSampleOverlap, final List<GenomeLoc> coverageIntervals) {
        super(dictionary, paddingFraction, minSampleOverlap);
        this.coverageIntervals = Utils.nonNull(coverageIntervals);
        genomicToBinMap = new TreeMap<>();
        for (int i = 0; i < coverageIntervals.size(); i++) {
            genomicToBinMap.put(coverageIntervals.get(i),i);
        }
        parser = new GenomeLocParser(this.dictionary);
    }

    public BinnedCNVDefragmenter(final SAMSequenceDictionary dictionary, final List<GenomeLoc> coverageIntervals) {
        this(dictionary, DEFAULT_PADDING_FRACTION, DEFAULT_SAMPLE_OVERLAP, coverageIntervals);
    }

    @Override
    protected SimpleInterval getPaddedRecordInterval(final SVCallRecord record) {
        Utils.nonNull(record);
        final GenomeLoc callStart = parser.createGenomeLoc(record.getContigA(), record.getPositionA(), record.getPositionA());
        final GenomeLoc callEnd = parser.createGenomeLoc(record.getContigA(), record.getPositionB(), record.getPositionB());

        //first interval that is equal to or "greater than" the call start, such that the start of the bin should match the call start, with a little wiggle room
        final Map.Entry<GenomeLoc, Integer> startBin = genomicToBinMap.ceilingEntry(callStart);
        if (startBin == null) {
            throw new UserException.BadInput("Call start " + callStart + " for  call " + record.getId() + " not found in model call intervals.");
        }
        final int callStartIndex = startBin.getValue();

        //last interval that is equal to or "less than" the call start, such that the end of the bin should match the call end
        final Map.Entry<GenomeLoc, Integer> endBin = genomicToBinMap.floorEntry(callEnd);
        if (endBin == null) {
            throw new UserException.BadInput("Call end " + callEnd + " for call " + record.getId() + " not found in model call intervals.");
        }
        final int callEndIndex = endBin.getValue();
        final int callBinLength = callEndIndex - callStartIndex + 1;
        if (callBinLength <= 0) {
            throw new UserException.BadInput("Copy number call at " + record.getContigA() + ":" + record.getPositionA() + "-"
                    + record.getPositionB() + " does not align with supplied model calling intervals. Use the filtered intervals input from GermlineCNVCaller for this cohort/model.");
        }

        final int paddedStartIndex = Math.max(callStartIndex - (int)Math.round(callBinLength * paddingFraction), 0);
        final int paddedCallStart;
        if (coverageIntervals.get(paddedStartIndex).getContig().equals(callStart.getContig())) {
            paddedCallStart = coverageIntervals.get(paddedStartIndex).getStart();
        } else {
            paddedCallStart = callStart.getStart();
        }

        final int paddedEndIndex = Math.min(callEndIndex + (int)Math.round(callBinLength * paddingFraction), genomicToBinMap.size() - 1);
        final int paddedCallEnd;
        if (coverageIntervals.get(paddedEndIndex).getContig().equals(callEnd.getContig())) {
            paddedCallEnd = coverageIntervals.get(paddedEndIndex).getEnd();
        } else {
            paddedCallEnd = callEnd.getEnd();
        }

        final int contigLength = dictionary.getSequence(record.getContigA()).getSequenceLength();
        return IntervalUtils.trimIntervalToContig(record.getContigA(), paddedCallStart, paddedCallEnd, contigLength);
    }

}
