package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.Set;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.COPY_NUMBER_FORMAT;

public class CNVDefragmenter extends SVClusterEngine<SVCallRecord> {

    protected static final double DEFAULT_SAMPLE_OVERLAP = 0.8;
    protected static final double DEFAULT_PADDING_FRACTION = 0.25;

    protected final double minSampleOverlap;
    protected final double paddingFraction;

    //for single-sample clustering case
    public CNVDefragmenter(final SAMSequenceDictionary dictionary, final double paddingFraction,
                           final double minSampleOverlap) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE, false,
                (new SVCollapser(SVCollapser.BreakpointSummaryStrategy.MIN_START_MAX_END))::collapse);
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
        if (!a.isDepthOnly() || !b.isDepthOnly()) return false;
        Utils.validate(a.getContigA().equals(a.getContigB()), "Variant A is depth-only but interchromosomal");
        Utils.validate(b.getContigA().equals(b.getContigB()), "Variant B is depth-only but interchromosomal");

        // Types match
        if (!a.getType().equals(b.getType())) return false;

        // Interval overlap
        if (!getPaddedRecordInterval(a.getContigA(), a.getPositionA(), a.getPositionB())
                .overlaps(getPaddedRecordInterval(b.getContigA(), b.getPositionA(), b.getPositionB()))) return false;

        // Sample overlap
        if (!hasSampleOverlap(a, b, minSampleOverlap)) {
            return false;
        }

        // In the single-sample case, match copy number strictly if we're looking at the same sample
        final Set<String> carriersA = getCopyNumberCarrierGenotypes(a).stream().map(Genotype::getSampleName).collect(Collectors.toSet());
        final Set<String> carriersB = getCopyNumberCarrierGenotypes(b).stream().map(Genotype::getSampleName).collect(Collectors.toSet());
        if (carriersA.size() == 1 && carriersB.size() == 1) {
            final String sampleA = carriersA.iterator().next();
            final String sampleB = carriersB.iterator().next();
            if (sampleA.equals(sampleB)) {
                final Genotype genotypeA = a.getGenotypes().get(sampleA);
                final Genotype genotypeB = b.getGenotypes().get(sampleB);
                if (genotypeA.hasExtendedAttribute(COPY_NUMBER_FORMAT) && genotypeB.hasExtendedAttribute(COPY_NUMBER_FORMAT)) {
                    final int copyNumberA = VariantContextGetters.getAttributeAsInt(genotypeA, COPY_NUMBER_FORMAT, 0);
                    final int copyNumberB = VariantContextGetters.getAttributeAsInt(genotypeB, COPY_NUMBER_FORMAT, 0);
                    if (copyNumberA != copyNumberB) {
                        return false;
                    }
                } else {
                    if (!genotypeA.getAlleles().equals(genotypeB.getAlleles())) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Max clusterable position depends on the other variant's length, which is unknown, so we calculate the size
     * of the largest possible event that would extend to the end of the contig.
     */
    @Override
    protected int getMaxClusterableStartingPosition(final SVCallRecord call) {
        final int contigLength = dictionary.getSequence(call.getContigA()).getSequenceLength();
        final int maxEventSize = (int) Math.ceil(contigLength - call.getPositionB());
        return Math.min((int) Math.ceil(call.getPositionB() + paddingFraction * (call.getLength() + maxEventSize)),
                dictionary.getSequence(call.getContigA()).getSequenceLength());
    }

    /**
     * Determine an overlap interval for clustering using padding specified at object construction
     * Returned interval represents the interval in which the start position of a new event must fall in order to be
     * added to the cluster (including the new event)
     * @return  an interval describing a cluster containing only this variant
     */
    protected SimpleInterval getPaddedRecordInterval(final String contig, final int start, final int end) {
        return new SimpleInterval(contig, start, end)
                .expandWithinContig((int) (paddingFraction * (end - start + 1)), dictionary);
    }

    public static double getDefaultPaddingFraction() {
        return DEFAULT_PADDING_FRACTION;
    }

    public static double getDefaultSampleOverlap() {
        return DEFAULT_SAMPLE_OVERLAP;
    }

    public double getPaddingFraction() {
        return paddingFraction;
    }
}
