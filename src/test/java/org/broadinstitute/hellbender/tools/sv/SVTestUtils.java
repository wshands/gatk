package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.LocatableClusterEngine;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine;
import org.broadinstitute.hellbender.tools.sv.cluster.SVCollapser;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.testng.Assert;

import java.util.*;
import java.util.function.Function;

public class SVTestUtils {
    public final static int chr1Length = 249250621;

    public final static SAMSequenceDictionary dict = new SAMSequenceDictionary(
            Arrays.asList(new SAMSequenceRecord("chr1", chr1Length),
                    new SAMSequenceRecord("chr10", 135534747),
                    new SAMSequenceRecord("chrX", 155270560)));

    private final static GenomeLocParser glParser = new GenomeLocParser(SVTestUtils.dict);
    public static final Function<Collection<SVCallRecord>, SVCallRecord> defaultCollapser =
            new SVCollapser(SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END)::collapse;

    public static SVClusterEngine<SVCallRecord> getNewDefaultSingleLinkageEngine() {
        final SVClusterEngine<SVCallRecord> engine = new SVClusterEngine<>(SVTestUtils.dict, LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, false, defaultCollapser);
        engine.setDepthOnlyParams(defaultDepthOnlyParameters);
        engine.setMixedParams(defaultMixedParameters);
        engine.setEvidenceParams(defaultEvidenceParameters);
        return engine;
    }

    public static SVClusterEngine<SVCallRecord> getNewDefaultMaxCliqueEngine() {
        final SVClusterEngine<SVCallRecord> engine = new SVClusterEngine<>(SVTestUtils.dict, LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE, false, defaultCollapser);
        engine.setDepthOnlyParams(defaultDepthOnlyParameters);
        engine.setMixedParams(defaultMixedParameters);
        engine.setEvidenceParams(defaultEvidenceParameters);
        return engine;
    }

    protected static final SVClusterEngine.ClusteringParameters defaultDepthOnlyParameters =
            new SVClusterEngine.DepthClusteringParameters(0.8, 0, 0);
    protected static final SVClusterEngine.ClusteringParameters defaultMixedParameters =
            new SVClusterEngine.MixedClusteringParameters(0.8, 1000, 0);
    protected static final SVClusterEngine.ClusteringParameters defaultEvidenceParameters =
            new SVClusterEngine.EvidenceClusteringParameters(0.5, 500, 0);

    public static final SVClusterEngine<SVCallRecord> defaultSingleLinkageEngine = getNewDefaultSingleLinkageEngine();
    public static final SVClusterEngine<SVCallRecord> defaultMaxCliqueEngine = getNewDefaultMaxCliqueEngine();

    public final static int start = 10001;

    public final static int length = 10000;

    public final static int length1b = (int)Math.round(defaultSingleLinkageEngine.getDepthOnlyParams().getReciprocalOverlap() * length);

    public final static int start1b = start + length - length1b;

    //separated from end of call1 by defragmenter padding (in bin space, according to bins defined by targetIntervals below)
    public final static int start2 = start + length*14/9;

    public final static int start3 = start + (int)Math.round(Math.floor((1- defaultSingleLinkageEngine.getDepthOnlyParams().getReciprocalOverlap())*length));

    //make intervals like xxxx----xxxx----xxxx----xxxx----xxxx
    public final static List<GenomeLoc> targetIntervals = new ArrayList<>(
            Arrays.asList(glParser.createGenomeLoc("chr1", 1, length/2),  //left edge call
                    glParser.createGenomeLoc("chr1", start, start + length/9),  //start of call1
                    glParser.createGenomeLoc("chr1", start + length*2/9, start + length*3/9),
                    glParser.createGenomeLoc("chr1", start + length*4/9, start + length*5/9),
                    glParser.createGenomeLoc("chr1", start + length*6/9, start + length*7/9),
                    glParser.createGenomeLoc("chr1", start + length*8/9, start + length),  //end of call1
                    glParser.createGenomeLoc("chr1", start + length*10/9, start + length*11/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start + length*12/9, start + length*13/9),  //padding for call1
                    glParser.createGenomeLoc("chr1", start2, start2 + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length*2/9, start2 + length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length*4/9, start2 + length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length*6/9, start2 + length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length*8/9, start2 + length),
                    glParser.createGenomeLoc("chr1", start2 + length, start2 + length + length/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*2/9, start2 + length+ length*3/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*4/9, start2 + length+ length*5/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*6/9, start2 + length+ length*7/9),
                    glParser.createGenomeLoc("chr1", start2 + length+ length*8/9, start2 + length+ length),
                    glParser.createGenomeLoc("chr1", chr1Length - 99, chr1Length),
                    glParser.createGenomeLoc("chr10", start, start + length - 1)))
            ;



    public final static GenotypeBuilder gb1 = new GenotypeBuilder("sample1", Collections.singletonList(Allele.create("<"+ GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)));
    public final static GenotypeBuilder gb2 = new GenotypeBuilder("sample1", Collections.singletonList(Allele.create("<"+ GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)));

    public final static Genotype sample1 = gb1.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 1).make();

    public final static Genotype sample1_CN0 = gb2.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0).make();

    public final static Map<String, Object> sample2AttributeMap = Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 3);

    public final static Map<String, Object> sample3AttributeMap = Collections.singletonMap(GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);

    public final static Genotype sample2 = GenotypeBuilder.create("sample2",
            Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_STRING_DUP +">", false)), sample2AttributeMap);

    public final static SVCallRecord rightEdgeCall = new SVCallRecord("rightEdgeCall", "chr1", chr1Length - 99, true,
            "chr1", chr1Length, true,
            StructuralVariantType.CNV, 100,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord leftEdgeCall = new SVCallRecord("leftEdgeCall", "chr1", 1, true,
            "chr1", length/2, true,
            StructuralVariantType.CNV, length/2,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord nonDepthOnly = new SVCallRecord("nonDepthOnly", "chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PE"),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1 = new SVCallRecord("call1", "chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1b = new SVCallRecord("call1b", "chr1", start1b, true,
            "chr1", start1b + length1b + 1, true,
            StructuralVariantType.CNV, length1b,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call2 = new SVCallRecord("call2", "chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call3 = new SVCallRecord("call3", "chr1", start2+length, true,
            "chr1", start2 + length + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call4_chr10 = new SVCallRecord("call4_chr10", "chr10", start, true,
            "chr10", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));

    public final static SVCallRecord call1_CN1 = new SVCallRecord("call1_CN1", "chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1));

    public final static SVCallRecord call2_CN0 = new SVCallRecord("call2_CN0", "chr1", start2, true,
            "chr1", start2 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1_CN0));

    public final static Genotype sample3 = GenotypeBuilder.create("sample3",
            Collections.singletonList(Allele.create("<"+GATKSVVCFConstants.SYMB_ALT_STRING_DEL +">", false)), sample3AttributeMap);

    public final static SVCallRecord sameBoundsSampleMismatch = new SVCallRecord("sameBoundsSampleMismatch", "chr1", start, true,
            "chr1", start + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
            Arrays.asList(sample1, sample3));

    public final static List<Genotype> threeGenotypes = Arrays.asList(sample1, sample2, sample3);

    public final static SVCallRecord inversion = new SVCallRecord("inversion", "chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.INV, 10001, Arrays.asList("SR", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthOnly = new SVCallRecord("depthOnly", "chr1", 10000, true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001, Collections.singletonList("depth"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL), Collections.singletonList(sample1));

    public static final SVCallRecord depthAndStuff = new SVCallRecord("depthAndStuff", "chr1", 10000, true, "chr1", 20000, true,
                    StructuralVariantType.CNV, 10001, Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff2 = new SVCallRecord("depthAndStuff2", "chr1", 10000 - defaultSingleLinkageEngine.getEvidenceParams().getWindow(), true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001 + defaultSingleLinkageEngine.getEvidenceParams().getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public static final SVCallRecord depthAndStuff3 = new SVCallRecord("depthAndStuff3", "chr1", 10000 + defaultSingleLinkageEngine.getEvidenceParams().getWindow(), true, "chr1", 20000, true,
            StructuralVariantType.CNV, 10001 - defaultSingleLinkageEngine.getEvidenceParams().getWindow(), Arrays.asList("depth", "PE"), Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP), Collections.singletonList(sample2));

    public final static SVCallRecord overlapsCall1 = new SVCallRecord("overlapsCall1", "chr1", start3, true,
            "chr1", start3 + length -1, true,
            StructuralVariantType.CNV, length,
            Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
            Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
            Arrays.asList(sample1, sample2));


    public static void assertEqualsExceptMembership(final SVCallRecord one, final SVCallRecord two) {
        if (one == two) return;
        Assert.assertEquals(one.getAlgorithms(), two.getAlgorithms());
        Assert.assertEquals(one.getAlleles(), two.getAlleles());

        // Exclude MEMBERS field
        final Map<String, Object> attributesOne = new HashMap<>(one.getAttributes());
        attributesOne.remove(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        final Map<String, Object> attributesTwo = new HashMap<>(two.getAttributes());
        attributesTwo.remove(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        Assert.assertEquals(attributesOne, attributesTwo);

        Assert.assertEquals(one.getPositionAInterval(), two.getPositionAInterval());
        Assert.assertEquals(one.getPositionBInterval(), two.getPositionBInterval());
        Assert.assertEquals(one.getStrandA(), two.getStrandA());
        Assert.assertEquals(one.getStrandB(), two.getStrandB());
        Assert.assertEquals(one.getId(), two.getId());
        Assert.assertEquals(one.getLength(), two.getLength());
        Assert.assertEquals(one.getType(), two.getType());
        Assert.assertEquals(one.getGenotypes().size(), two.getGenotypes().size());
        for (int i = 0; i < one.getGenotypes().size(); i++) {
            VariantContextTestUtils.assertGenotypesAreEqual(one.getGenotypes().get(i), two.getGenotypes().get(i));
        }
    }

    public static void assertContainsAll(final GenotypesContext one, final GenotypesContext two) {
        Assert.assertTrue(one.size() >= two.size());
        if ( !one.isEmpty()) {
            Assert.assertTrue(one.getSampleNames().containsAll(two.getSampleNames()), "sample names set");
            final Set<String> samples = two.getSampleNames();
            for ( final String sample : samples ) {
                VariantContextTestUtils.assertGenotypesAreEqual(one.get(sample), two.get(sample));
            }
        }
    }
}
