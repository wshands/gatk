package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.function.Supplier;

public class CollatingInterval implements Feature, Comparable<CollatingInterval> {
    private final SAMSequenceRecord contig;
    private final int start;
    private final int end;

    public CollatingInterval( final SAMSequenceDictionary dict, final Locatable loc ) {
        contig = dict.getSequence(loc.getContig());
        start = loc.getStart();
        end = loc.getEnd();
        validate(() -> loc.getContig() + " not in dictionary");
    }

    public CollatingInterval( final SAMSequenceDictionary dict, final String contigName,
                              final int start, final int end ) {
        this.contig = dict.getSequence(contigName);
        this.start = start;
        this.end = end;
        validate(() -> contigName + " not in dictionary");
    }

    public CollatingInterval( final SAMSequenceRecord contig, final int start, final int end ) {
        this.contig = contig;
        this.start = start;
        this.end = end;
        validate(() -> "null contig");
    }

    public CollatingInterval( final SAMSequenceDictionary dict, final DataInputStream dis )
        throws IOException {
        final int contigId = dis.readInt();
        contig = dict.getSequence(contigId);
        start = dis.readInt();
        end = dis.readInt();
        validate(() -> "dictionary does not contain a contig with id " + contigId);
    }

    @Override public String getContig() { return contig.getSequenceName(); }
    @Override public int getStart() { return start; }
    @Override public int getEnd() { return end; }
    @Override public boolean overlaps( final Locatable that ) {
        return contigsMatch(that) && start <= that.getEnd() && that.getStart() <= end;
    }
    @Override public boolean contains( final Locatable that ) {
        return contigsMatch(that) && that.getStart() >= start && that.getEnd() <= end;
    }
    @Override public boolean contigsMatch( final Locatable that ) {
        final String thisContig = contig.getSequenceName();
        final String thatContig = that.getContig();
        return thisContig == thatContig || thisContig.equals(thatContig);
    }

    /**
     * Compares contig, start, and stop, in that order.
     * The contig comparison is done relative to some dictionary's order.
     */
    @Override public int compareTo( final CollatingInterval that ) {
        int result = Integer.compare(this.contig.getSequenceIndex(), that.contig.getSequenceIndex());
        if ( result == 0 ) {
            result = Integer.compare(this.start, that.start);
            if ( result == 0 ) {
                result = Integer.compare(this.end, that.end);
            }
        }
        return result;
    }

    @Override public boolean equals( final Object obj ) {
        if ( this == obj ) return true;
        if ( !(obj instanceof CollatingInterval) ) return false;
        final CollatingInterval that = (CollatingInterval)obj;
        return contigsMatch(that) && this.start == that.start && this.end == that.end;
    }

    @Override public int hashCode() {
        return 241*(241*(241*contig.getSequenceIndex() + start) + end);
    }

    @Override public String toString() {
        return contig.getSequenceName() + ":" + start + "-" + end;
    }

    /**
     * Returns a value less than 0 for a locus earlier in the genome than any in this interval,
     * the value 0 for a locus that overlaps this interval,
     * and a value greater than 0 for a locus later than any any in this interval.
     */
    public int compareLocus( final SAMSequenceRecord thatContig, final int thatPosition ) {
        int result = Integer.compare(thatContig.getSequenceIndex(), contig.getSequenceIndex());
        if ( result != 0 ) {
            return result;
        }
        if ( thatPosition < start ) {
            return -1;
        }
        if ( thatPosition > end ) {
            return 1;
        }
        return 0;
    }

    /**
     * Upstream means not overlapping, and
     * 1) this is on an earlier (relative to some dictionary's order) contig than that, or
     * 2) if on the same contig, this ends earlier than that starts
     */
    public boolean isUpstreamOf( final CollatingInterval that ) {
        final int thisContigId = contig.getSequenceIndex();
        final int thatContigId = that.contig.getSequenceIndex();
        if ( thisContigId < thatContigId ) {
            return true;
        }
        if ( thisContigId == thatContigId ) {
            return end < that.getStart();
        }
        return false;
    }

    /**
     * If the contigs are not the same, the laterEnding interval is the one with the later (relative
     * to some dictionary) contig.
     * If the contigs are the same, the laterEnding interval is the one with the greater end position.
     */
    public static CollatingInterval laterEnding( final CollatingInterval interval1,
                                                 final CollatingInterval interval2 ) {
        final int contigId1 = interval1.contig.getSequenceIndex();
        final int contigId2 = interval2.contig.getSequenceIndex();
        if ( contigId1 == contigId2 ) {
            return interval1.end > interval2.end ? interval1 : interval2;
        }
        if ( contigId1 > contigId2 ) {
            return interval1;
        }
        return interval2;
    }

    public void write( final DataOutputStream dos ) throws IOException {
        dos.writeInt(contig.getSequenceIndex());
        dos.writeInt(start);
        dos.writeInt(end);
    }

    private void validate( final Supplier<String> badContigMessage ) {
        if ( contig == null ) {
            throw new GATKException(badContigMessage.get());
        }
        final int sequenceLength = contig.getSequenceLength();
        if ( start < 1 || start > sequenceLength ) {
            throw new GATKException("starting coordinate " + start + " is not within contig bounds");
        }
        if ( end < start || end > sequenceLength ) {
            throw new GATKException("ending coordinate " + end +
                    " is less than start or greater than contig length");
        }
    }
}
