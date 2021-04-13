package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVSegmentVariantComposer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class SVCallRecordUtils {

    /**
     * Create a variant from a call for VCF interoperability
     *
     * @param call variant to convert
     * @return
     */
    public static VariantContextBuilder getVariantBuilder(final SVCallRecord call) {
        Utils.nonNull(call);
        final int end;
        if (call.getType().equals(StructuralVariantType.INS) || call.getType().equals(StructuralVariantType.BND)) {
            end = call.getPositionA();
        } else {
            end = call.getPositionB();
        }
        final int end2;
        if (call.getType().equals(StructuralVariantType.INS)) {
            end2 = call.getPositionA();
        } else {
            end2 = call.getPositionB();
        }

        final List<Allele> altAlleles = call.getAlleles().stream()
                .filter(a -> !a.isNoCall() && !a.isReference())
                .distinct()
                .collect(Collectors.toList());
        final List<Allele> refAlleles = call.getAlleles().stream()
                .filter(a -> !a.isNoCall() && a.isReference())
                .distinct()
                .collect(Collectors.toList());
        Utils.validate(refAlleles.size() <= 1, "Encountered multiple reference alleles");
        final int numAlleles = refAlleles.size() + altAlleles.size();
        final List<Allele> alleles = new ArrayList<>(numAlleles);
        if (refAlleles.isEmpty()) {
            // If no reference allele, use N
            alleles.add(Allele.REF_N);
        } else {
            alleles.add(refAlleles.get(0));
        }
        alleles.addAll(altAlleles);

        final VariantContextBuilder builder = new VariantContextBuilder(call.getId(), call.getContigA(), call.getPositionA(),
                end, alleles);
        builder.id(call.getId());
        builder.attributes(call.getAttributes());
        builder.attribute(VCFConstants.END_KEY, end);
        builder.attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, call.getContigB());
        builder.attribute(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        builder.attribute(GATKSVVCFConstants.SVLEN, call.getLength());
        builder.attribute(GATKSVVCFConstants.SVTYPE, call.getType());
        builder.attribute(GATKSVVCFConstants.STRANDS_ATTRIBUTE, getStrandString(call));
        builder.attribute(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, call.getAlgorithms());
        builder.genotypes(call.getGenotypes());
        return builder;
    }

    public static GenotypesContext fillMissingSamplesWithGenotypes(final GenotypesContext genotypes,
                                                                   final List<Allele> alleles,
                                                                   final Set<String> samples,
                                                                   final Map<String,Object> attributes) {
        Utils.nonNull(genotypes);
        Utils.nonNull(samples);
        final Set<String> missingSamples = Sets.difference(samples, genotypes.getSampleNames());
        if (missingSamples.isEmpty()) {
            return genotypes;
        }
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size() + missingSamples.size());
        newGenotypes.addAll(genotypes);
        for (final String sample : missingSamples) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            if (attributes != null) {
                genotypeBuilder.attributes(attributes);
            }
            genotypeBuilder.alleles(alleles);
            newGenotypes.add(genotypeBuilder.make());
        }
        return GenotypesContext.copy(newGenotypes);
    }

    public static SVCallRecord copyCallWithNewGenotypes(final SVCallRecord record, final GenotypesContext genotypes) {
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(), record.getAlleles(),
                genotypes, record.getAttributes());
    }

    public static SVCallRecordWithEvidence copyCallWithNewGenotypes(final SVCallRecordWithEvidence record, final GenotypesContext genotypes) {
        return new SVCallRecordWithEvidence(record.getId(), record.getContigA(), record.getPositionA(), record.getStrandA(), record.getContigB(),
                record.getPositionB(), record.getStrandB(), record.getType(), record.getLength(), record.getAlgorithms(), record.getAlleles(),
                genotypes, record.getAttributes(), record.getStartSplitReadSites(), record.getEndSplitReadSites(), record.getDiscordantPairs(),
                record.getCopyNumberDistribution());
    }

    public static VariantContextBuilder createBuilderWithEvidence(final SVCallRecordWithEvidence call, boolean updateSplitRead, boolean updateDiscordantPair) {
        final VariantContextBuilder builder = getVariantBuilder(call);
        final boolean includeEvidence = !call.isDepthOnly();
        final SplitReadSite startSplitReadCounts = includeEvidence ? getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getPositionA()) : null;
        final SplitReadSite endSplitReadCounts = includeEvidence ? getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getPositionB()) : null;
        final Map<String,Integer> discordantPairCounts = includeEvidence ? getDiscordantPairCountsMap(call.getDiscordantPairs()) : null;
        final List<Genotype> genotypes = builder.getGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(genotypes.size());
        for (final Genotype genotype : genotypes) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            if (includeEvidence) {
                if (updateSplitRead) {
                    genotypeBuilder.attribute(GATKSVVCFConstants.START_SPLIT_READ_COUNT_ATTRIBUTE, startSplitReadCounts.getCount(sample));
                    genotypeBuilder.attribute(GATKSVVCFConstants.END_SPLIT_READ_COUNT_ATTRIBUTE, endSplitReadCounts.getCount(sample));
                }
                if (updateDiscordantPair) {
                    genotypeBuilder.attribute(GATKSVVCFConstants.DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
                }
            }
            newGenotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(newGenotypes);
        return builder;
    }

    private static String getStrandString(final SVCallRecord call) {
        return getStrandString(call.getStrandA()) + getStrandString(call.getStrandB());
    }

    private static String getStrandString(final boolean strand) {
        return strand ? SVCallRecord.STRAND_PLUS : SVCallRecord.STRAND_MINUS;
    }

    private static SplitReadSite getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .findAny()
                .orElse(new SplitReadSite(pos, Collections.emptyMap()));
    }

    private static Map<String,Integer> getDiscordantPairCountsMap(final Collection<DiscordantPairEvidence> discordantPairs) {
        Utils.nonNull(discordantPairs);
        return discordantPairs.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
    }

    public static <T extends SVLocatable> Comparator<T> getSVLocatableComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareSVLocatables(o1, o2, dictionary);
    }

    public static <T extends SVLocatable> Comparator<T> getSVLocatableComparatorByEnds(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareSVLocatablesByEnds(o1, o2, dictionary);
    }

    public static <T extends SVCallRecord> Comparator<T> getCallComparator(final SAMSequenceDictionary dictionary) {
        return (o1, o2) -> compareCalls(o1, o2, dictionary);
    }

    public static int compareSVLocatablesByEnds(final SVLocatable first, final SVLocatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        // First locus
        final Comparator<Locatable> locatableComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        // Second locus
        final int compareB = locatableComparator.compare(new SimpleInterval(first.getContigB(), first.getPositionB(), first.getPositionB()),
                new SimpleInterval(second.getContigB(), second.getPositionB(), second.getPositionB()));
        return compareB;
    }

    public static int compareSVLocatables(final SVLocatable first, final SVLocatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        // First locus
        final Comparator<Locatable> locatableComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final int compareA = locatableComparator.compare(new SimpleInterval(first.getContigA(), first.getPositionA(), first.getPositionA()),
                new SimpleInterval(second.getContigA(), second.getPositionA(), second.getPositionA()));
        if (compareA != 0) return compareA;
        // Second locus
        final int compareB = locatableComparator.compare(new SimpleInterval(first.getContigB(), first.getPositionB(), first.getPositionB()),
                new SimpleInterval(second.getContigB(), second.getPositionB(), second.getPositionB()));
        return compareB;
    }

    public static int compareCalls(final SVCallRecord first, final SVCallRecord second, final SAMSequenceDictionary dictionary) {
        final int compareLocatables = compareSVLocatables(first, second, dictionary);
        if (compareLocatables != 0) return compareLocatables;

        //Strands
        final int compareStartStrand = Boolean.compare(first.getStrandA(), second.getStrandA());
        if (compareStartStrand != 0) return compareStartStrand;
        final int compareEndStrand = Boolean.compare(first.getStrandB(), second.getStrandB());
        if (compareEndStrand != 0) return compareEndStrand;

        // Length
        final int compareLength = Integer.compare(first.getLength(), second.getLength());
        if (compareLength != 0) return compareLength;

        // Type
        final int compareType = first.getType().compareTo(second.getType());
        return compareType;
    }

    public static boolean isValidSize(final SVCallRecord call, final int minEventSize) {
        // Treat inter-chromosomal events as always valid
        if (!call.getContigA().equals(call.getContigB())) {
            return true;
        }
        // BNDs have undefined nominal length, so use coordinates
        if (call.getType().equals(StructuralVariantType.BND)) {
            return Math.abs(call.getPositionB() - call.getPositionA()) >= minEventSize;
        }
        // Otherwise use length
        return call.getLength() >= minEventSize;
    }

    public static <T> boolean intervalIsIncluded(final SVCallRecord call, final Map<String, IntervalTree<T>> includedIntervalTreeMap,
                                                 final double minOverlapFraction, final boolean requireBreakendOverlap) {
        final IntervalTree<T> startTree = includedIntervalTreeMap.get(call.getContigA());
        // Contig A included
        if (startTree == null) {
            return false;
        }
        final IntervalTree<T> endTree = includedIntervalTreeMap.get(call.getContigB());
        // Contig B included
        if (endTree == null) {
            return false;
        }
        // Breakends both included, if required
        if (requireBreakendOverlap && !(startTree.overlappers(call.getPositionA(), call.getPositionA() + 1).hasNext()
                && endTree.overlappers(call.getPositionB(), call.getPositionB() + 1).hasNext())) {
            return false;
        }
        // Can include all interchromosomal variants at this point
        if (!call.getContigA().equals(call.getContigB())) {
            return true;
        }
        // Check overlap fraction
        final long overlapLength = totalOverlap(call.getPositionA(), call.getPositionB(), startTree);
        final double overlapFraction = overlapLength / (double) (call.getPositionB() - call.getPositionA());
        return overlapFraction >= minOverlapFraction;
    }

    private static <T> long totalOverlap(final int start, final int end, final IntervalTree<T> tree) {
        final Iterator<IntervalTree.Node<T>> iter = tree.overlappers(start, end);
        long overlap = 0;
        while (iter.hasNext()) {
            final IntervalTree.Node<T> node = iter.next();
            overlap += intersectionLength(start, end, node.getStart(), node.getEnd());
        }
        return overlap;
    }

    private static long intersectionLength(final int start1, final int end1, final int start2, final int end2) {
        return Math.max(0, Math.min(end1, end2) - Math.max(start1, start2) + 1);
    }

    public static Stream<SVCallRecord> convertInversionsToBreakends(final SVCallRecord call) {
        if (!call.getType().equals(StructuralVariantType.INV)) {
            return Stream.of(call);
        }
        Utils.validateArg(isIntrachromosomal(call), "Inversion is not intrachromosomal");
        final SVCallRecord positiveBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), true, call.getContigB(), call.getPositionB(), true, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getAlleles(), call.getGenotypes(), call.getAttributes());
        final SVCallRecord negativeBreakend = new SVCallRecord(call.getId(), call.getContigA(),
                call.getPositionA(), false, call.getContigB(), call.getPositionB(), false, StructuralVariantType.BND, -1,
                call.getAlgorithms(), call.getAlleles(), call.getGenotypes(), call.getAttributes());
        return Stream.of(positiveBreakend, negativeBreakend);
    }

    public static void validateCoordinates(final SVCallRecord call) {
        Utils.nonNull(call);
        Utils.validateArg(call.getPositionA() >= 1, "Call start non-positive");
        if (isIntrachromosomal(call)) {
            Utils.validateArg(call.getPositionA() <= call.getPositionB(), "Second position before end on same contig");
        } else {
            Utils.validateArg(call.getPositionB() >= 1, "Call second position non-positive");
        }
    }

    public static void validateCoordinatesWithDictionary(final SVCallRecord call, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);
        validateCoordinates(call);
        final SAMSequenceRecord contigARecord = dictionary.getSequence(call.getContigA());
        Utils.validateArg(contigARecord != null, "Call first contig " + call.getContigA() + " not in dictionary");
        final SAMSequenceRecord contigBRecord = dictionary.getSequence(call.getContigB());
        Utils.validateArg(contigBRecord != null, "Call second contig " + call.getContigB() + " not in dictionary");
        Utils.validateArg(call.getPositionA() <= contigARecord.getSequenceLength(), "Call first position greater than contig length");
        Utils.validateArg(call.getPositionB() <= contigBRecord.getSequenceLength(), "Call second position greater than contig length");
    }

    public static boolean isIntrachromosomal(final SVCallRecord call) {
        return call.getContigA().equals(call.getContigB());
    }

    private static SplitReadSite createSplitReadSite(final VariantContext variant, final int position, final String splitReadCountKey) {
        final Map<String,Integer> counts = variant.getGenotypes().stream()
                .collect(Collectors.toMap(g -> g.getSampleName(), g -> (Integer) g.getExtendedAttribute(splitReadCountKey)));
        return new SplitReadSite(position, counts);
    }

    public static SVCallRecord create(final VariantContext variant) {
        return create(variant, true);
    }

    public static SVCallRecord create(final VariantContext variant, boolean keepVariantAttributes) {
        Utils.nonNull(variant);
        final String id = variant.getID();
        final String contigA = variant.getContig();
        final int positionA = variant.getStart();

        final StructuralVariantType type = inferStructuralVariantType(variant);
        final List<String> algorithms = getAlgorithms(variant);
        final String strands = getStrands(variant, type);
        final boolean strand1 = strands.startsWith(SVCallRecord.STRAND_PLUS);
        final boolean strand2 = strands.endsWith(SVCallRecord.STRAND_PLUS);
        final int length = getLength(variant);
        final Map<String, Object> attributes = keepVariantAttributes ? variant.getAttributes() : Collections.emptyMap();

        final String contigB;
        final int positionB;
        if (type.equals(StructuralVariantType.BND)) {
            // If END2 and CONTIG2 are both defined, use those.
            // If neither is defined, use start contig and position.
            // If only CONTIG2 is defined, END2 is taken as END
            // Having only END2 but not CONTIG2 is unacceptable
            final boolean hasContig2 = variant.hasAttribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE);
            final boolean hasEnd2 = variant.hasAttribute(GATKSVVCFConstants.END2_ATTRIBUTE);
            if (!(hasContig2 && hasEnd2)) {
                throw new UserException.BadInput("Attributes " + GATKSVVCFConstants.END2_ATTRIBUTE +
                        " and " + GATKSVVCFConstants.CONTIG2_ATTRIBUTE + " are required for BND records.");
            }
            contigB = variant.getAttributeAsString(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, null);
            positionB = variant.getAttributeAsInt(GATKSVVCFConstants.END2_ATTRIBUTE, 0);
        } else {
            contigB = contigA;
            positionB = variant.getEnd();
        }
        return new SVCallRecord(id, contigA, positionA, strand1, contigB, positionB, strand2, type, length, algorithms, variant.getAlleles(), variant.getGenotypes(), attributes);
    }

    public static int getLength(final VariantContext variant) {
        Utils.nonNull(variant);
        if (variant.hasAttribute(GATKSVVCFConstants.SVLEN)) {
            return variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        }
        return variant.getEnd() - variant.getStart(); // TODO +1?
    }

    public static List<String> getAlgorithms(final VariantContext variant) {
        Utils.nonNull(variant);
        final List<String> algorithmsAttr = variant.getAttributeAsStringList(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, null);
        if (algorithmsAttr != null) {
            return algorithmsAttr;
        }
        final StructuralVariantType type = variant.getStructuralVariantType();
        if (type.equals(StructuralVariantType.DEL) || type.equals(StructuralVariantType.DUP)) {
            return Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
        }
        throw new UserException.BadInput("Non-DEL/DUP variant missing " + GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE + " attribute.");
    }

    public static String getStrands(final VariantContext variant, final StructuralVariantType type) {
        Utils.nonNull(variant);
        Utils.nonNull(type);
        final String strandsAttr = variant.getAttributeAsString(GATKSVVCFConstants.STRANDS_ATTRIBUTE, null);
        if (strandsAttr == null) {
            if (type.equals(StructuralVariantType.DEL) || type.equals(StructuralVariantType.INS)) {
                return SVCallRecord.STRAND_PLUS + SVCallRecord.STRAND_MINUS;
            } else if (type.equals(StructuralVariantType.DUP)) {
                return SVCallRecord.STRAND_MINUS + SVCallRecord.STRAND_PLUS;
            } else {
                throw new UserException.BadInput("Record of type " + type.name() + " missing required attribute "
                        + GATKSVVCFConstants.STRANDS_ATTRIBUTE);
            }
        }
        if (strandsAttr.length() != 2) {
            throw new IllegalArgumentException("Strands field is not 2 characters long");
        }
        final String strand1Char = strandsAttr.substring(0, 1);
        if (!strand1Char.equals(SVCallRecord.STRAND_PLUS) && !strand1Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final String strand2Char = strandsAttr.substring(1, 2);
        if (!strand2Char.equals(SVCallRecord.STRAND_PLUS) && !strand2Char.equals(SVCallRecord.STRAND_MINUS)) {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        return strandsAttr;
    }

    public static StructuralVariantType inferStructuralVariantType(final VariantContext variant) {
        final StructuralVariantType type = variant.getStructuralVariantType();
        if (type != null) {
            return type;
        }
        final List<Allele> alleles = variant.getAlternateAlleles();
        Utils.validate(!alleles.isEmpty(), "Variant missing alt allele");
        if (alleles.size() == 2 && alleles.contains(GATKSVVCFConstants.DEL_ALLELE) && alleles.contains(GATKSVVCFConstants.DUP_ALLELE)) {
            return StructuralVariantType.CNV;
        }
        Utils.validate(alleles.size() == 1, "Non-CNV multiallelic variants not supported");
        final Allele allele = alleles.get(0);
        Utils.validate(allele.isSymbolic(), "Expected symbolic alt allele");
        return StructuralVariantType.valueOf(allele.getDisplayString().replace("<", "").replace(">", ""));
    }

    /**
     *
     * @param variant single-sample variant from a gCNV segments VCF
     * @param minQuality drop events with quality lower than this
     * @return
     */
    public static SVCallRecord createDepthOnlyFromGCNVWithOriginalGenotypes(final VariantContext variant, final double minQuality) {
        Utils.nonNull(variant);
        if (variant.getGenotypes().size() == 1) {
            //only cluster good variants
            final Genotype g = variant.getGenotypes().get(0);
            if (g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                    || VariantContextGetters.getAttributeAsInt(g, GermlineCNVSegmentVariantComposer.QS, 0) < minQuality
                    || isNullCall(g)) {
                return null;
            }
        }

        final SVCallRecord baseRecord;
        if (variant.getReference().equals(Allele.REF_N)) { //old segments VCFs had ref Ns and genotypes that didn't reflect ploidy accurately
            baseRecord = createVariantFromLegacyGCNV(variant);
        } else {
            baseRecord = create(variant, true);
        }
        final List<Genotype> nonRefGenotypes = baseRecord.getGenotypes().stream()
                .filter(g -> !(g.isHomRef() || (g.isNoCall() && !g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))))
                .collect(Collectors.toList());
        return copyCallWithNewGenotypes(baseRecord, GenotypesContext.copy(nonRefGenotypes));
    }

    private static SVCallRecord createVariantFromLegacyGCNV(final VariantContext variant) {
        boolean isDel = false;
        for (final Genotype g : variant.getGenotypes()) {
            if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DEL_ALLELE))) {
                isDel = true;
            } else if (g.getAlleles().stream().anyMatch(a -> a.equals(GATKSVVCFConstants.DUP_ALLELE))) {
                isDel = false;
            } else if (g.getAlleles().stream().allMatch(a -> a.isNoCall())) {
                if (g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
                    isDel = (Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()) < g.getPloidy());
                } else {
                    throw new IllegalStateException("Genotype for sample " + g.getSampleName() + " at " + variant.getContig() + ":" + variant.getStart() + " had no CN attribute and will be dropped.");
                }
            } else {
                throw new IllegalArgumentException("Segment VCF schema expects <DEL>, <DUP>, and no-call allele, but found "
                        + g.getAllele(0) + " at " + variant.getContig() + ":" + variant.getStart());
            }
        }
        final boolean startStrand = isDel ? true : false;
        final boolean endStrand = isDel ? false : true;
        final StructuralVariantType type;
        final List<Allele> alleles;
        if (!variant.getReference().equals(Allele.REF_N) && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DUP_ALLELE)
                && variant.getAlternateAlleles().contains(GATKSVVCFConstants.DEL_ALLELE)) {
            type = StructuralVariantType.CNV;
            alleles = variant.getAlleles();
        } else if (isDel) {
            type = StructuralVariantType.DEL;
            alleles = Lists.newArrayList(variant.getReference(), GATKSVVCFConstants.DEL_ALLELE);
        } else {
            type = StructuralVariantType.DUP;
            alleles = Lists.newArrayList(variant.getReference(), GATKSVVCFConstants.DUP_ALLELE);
        }

        final String id = variant.getID();
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final int end = variant.getEnd();
        final int length = end - start;
        final List<String> algorithms = Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM);
        return new SVCallRecord(id, startContig, start, startStrand, startContig, end, endStrand, type, length, algorithms, alleles, variant.getGenotypes(), variant.getAttributes());
    }

    /**
     *
     * @param g
     * @return true if this is a call on a missing contig
     */
    private static boolean isNullCall(final Genotype g) {
        return g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)
                && VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0) == 0
                && g.isNoCall();

    }

    public static boolean isAltGenotype(final Genotype g) {
        return g.getAlleles().stream().anyMatch(a -> !a.isNoCall() && !a.isReference());
    }
}
