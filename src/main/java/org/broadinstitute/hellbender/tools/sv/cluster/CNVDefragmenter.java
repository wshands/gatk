package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedHashSet;
import java.util.Set;

public class CNVDefragmenter extends SVClusterEngine<SVCallRecord> {

    protected static final double DEFAULT_SAMPLE_OVERLAP = 0.8;
    protected static final double DEFAULT_PADDING_FRACTION = 0.25;

    protected final double minSampleOverlap;
    protected final double paddingFraction;

    //for single-sample clustering case
    public CNVDefragmenter(final SAMSequenceDictionary dictionary, final double paddingFraction,
                           final double minSampleOverlap) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE, false,
                (new CNVCollapser(SVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END))::collapse);
        this.minSampleOverlap = minSampleOverlap;
        this.paddingFraction = paddingFraction;
    }

    public CNVDefragmenter(final SAMSequenceDictionary dictionary) {
        this(dictionary, DEFAULT_PADDING_FRACTION, DEFAULT_SAMPLE_OVERLAP);
    }

    /**
     * Determine if two variants should cluster based on their padded intervals and genotyped samples
     * @param a
     * @param b
     * @return true if the two variants should be in the same cluster
     */
    @Override
    protected boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        // Only do clustering on depth-only variants
        if (!isDepthOnlyCall(a) || !isDepthOnlyCall(b)) return false;
        Utils.validate(a.getContigA().equals(a.getContigB()), "Variant A is depth-only but interchromosomal");
        Utils.validate(b.getContigA().equals(b.getContigB()), "Variant B is depth-only but interchromosomal");

        // Types match
        if (!a.getType().equals(b.getType())) return false;

        // Interval overlap
        if (!getPaddedRecordInterval(a).overlaps(getPaddedRecordInterval(b))) return false;

        // Sample overlap
        final Set<String> sharedSamples = new LinkedHashSet<>(a.getCalledSamples());
        sharedSamples.retainAll(b.getCalledSamples());
        final double sampleOverlap = Math.min(sharedSamples.size() / (double) a.getCalledSamples().size(), sharedSamples.size() / (double) b.getCalledSamples().size());
        if (sampleOverlap < minSampleOverlap) return false;

        // In the single-sample case, match copy number strictly if we're looking at the same sample
        if (a.getGenotypes().size() == 1 && b.getGenotypes().size() == 1
                && a.getGenotypes().get(0).getSampleName().equals(b.getGenotypes().get(0).getSampleName())
                && a.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)
                && b.getGenotypes().get(0).hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) &&
                !(a.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).equals(b.getGenotypes().get(0).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)))) {
            return false;
        }
        return true;
    }

    /**
     * Determine an overlap interval for clustering using padding specified at object construction
     * Returned interval represents the interval in which the start position of a new event must fall in order to be
     * added to the cluster (including the new event)
     * @param record  new event to be clustered
     * @return  an interval describing a cluster containing only this variant
     */
    protected SimpleInterval getPaddedRecordInterval(final SVCallRecord record) {
        return new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB())
                .expandWithinContig((int) (paddingFraction * record.getLength()), dictionary);
    }

    public static double getDefaultPaddingFraction() {
        return DEFAULT_PADDING_FRACTION;
    }

    public static double getDefaultSampleOverlap() {
        return DEFAULT_SAMPLE_OVERLAP;
    }

    public static boolean isDepthOnlyCall(final SVCallRecord record) {
        if (record.getAlgorithms().isEmpty()) return false;
        for (final String alg : record.getAlgorithms()) {
            if (!alg.equals(GATKSVVCFConstants.DEPTH_ALGORITHM)) return false;
        }
        return true;
    }

    public double getPaddingFraction() {
        return paddingFraction;
    }
}
