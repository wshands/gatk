package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.COPY_NUMBER_FORMAT;

public class CNVDefragmenter extends SVClusterEngine<SVCallRecord> {

    public static final double DEFAULT_SAMPLE_OVERLAP = 0.8;
    public static final double DEFAULT_PADDING_FRACTION = 0.25;

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

    /**
     * Determine if two variants should cluster based on their padded intervals and genotyped samples
     * @param a
     * @param b
     * @return true if the two variants should be in the same cluster
     */
    @Override
    boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        // Only do clustering on depth-only CNVs
        if (!a.isDepthOnly() || !b.isDepthOnly()) return false;
        if (!a.isCNV() || !b.isCNV()) return false;
        Utils.validate(a.getContigA().equals(a.getContigB()), "Variant A is a CNV but interchromosomal");
        Utils.validate(b.getContigA().equals(b.getContigB()), "Variant B is a CNV but interchromosomal");

        // Types match
        if (a.getType() != b.getType()) return false;

        // Interval overlap
        if (!getPaddedRecordInterval(a.getContigA(), a.getPositionA(), a.getPositionB())
                .overlaps(getPaddedRecordInterval(b.getContigA(), b.getPositionA(), b.getPositionB()))) return false;

        // Sample overlap
        if (!hasSampleOverlap(a, b, minSampleOverlap)) {
            return false;
        }

        // In the single-sample case, match copy number strictly if we're looking at the same sample
        // TODO repeated check for CN attributes in hasSampleOverlap and getBestAvailableCarrierSamples
        final List<String> carriersA = getBestAvailableCarrierSamples(a);
        final List<String> carriersB = getBestAvailableCarrierSamples(b);
        if (carriersA.size() == 1 && carriersA.equals(carriersB)) {
            final Genotype genotypeA = a.getGenotypes().get(carriersA.get(0));
            final Genotype genotypeB = b.getGenotypes().get(carriersB.get(0));
            if (genotypeA.hasExtendedAttribute(COPY_NUMBER_FORMAT) && genotypeB.hasExtendedAttribute(COPY_NUMBER_FORMAT)) {
                final int copyNumberA = VariantContextGetters.getAttributeAsInt(genotypeA, COPY_NUMBER_FORMAT, 0);
                final int copyNumberB = VariantContextGetters.getAttributeAsInt(genotypeB, COPY_NUMBER_FORMAT, 0);
                final int copyNumberDeltaA = genotypeA.getPloidy() - copyNumberA;
                final int copyNumberDeltaB = genotypeB.getPloidy() - copyNumberB;
                if (copyNumberDeltaA != copyNumberDeltaB) {
                    return false;
                }
            } else {
                final List<Allele> sortedAllelesA = genotypeA.getAlleles().stream().sorted().collect(Collectors.toList());
                final List<Allele> sortedAllelesB = genotypeB.getAlleles().stream().sorted().collect(Collectors.toList());
                if (!sortedAllelesA.equals(sortedAllelesB)) {
                    return false;
                }
            }
        }
        return true;
    }

    private List<String> getBestAvailableCarrierSamples(final SVCallRecord record) {
        if (hasDefinedCopyNumbers(record.getGenotypes())) {
            return getCopyNumberCarrierGenotypes(record).stream().map(Genotype::getSampleName).collect(Collectors.toList());
        } else {
            return getGenotypedCarrierGenotypes(record, new HashSet<>(record.getAltAlleles())).stream().map(Genotype::getSampleName).collect(Collectors.toList());
        }
    }

    /**
     * Max clusterable position depends on the other variant's length, which is unknown, so we calculate the size
     * of the largest possible event that would extend to the end of the contig.
     */
    @Override
    protected int getMaxClusterableStartingPosition(final SVCallRecord call) {
        final int contigLength = dictionary.getSequence(call.getContigA()).getSequenceLength();
        final int maxTheoreticalStart = (int) Math.floor((call.getPositionB() + paddingFraction * (call.getLength() + contigLength)) / (1.0 + paddingFraction));
        return Math.min(maxTheoreticalStart, dictionary.getSequence(call.getContigA()).getSequenceLength());
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
}
