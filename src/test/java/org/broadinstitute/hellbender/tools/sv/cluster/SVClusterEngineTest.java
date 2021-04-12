package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.sv.cluster.LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
import static org.broadinstitute.hellbender.tools.sv.cluster.LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE;

public class SVClusterEngineTest {

    private final SVClusterEngine<SVCallRecord> engine = SVTestUtils.defaultSingleLinkageEngine;

    @BeforeTest
    public void initializeClusterEngine() {
        engine.add(SVTestUtils.call1);
    }

    @Test
    public void testCollapser() {
        //depth only and depthAndStuff have same bounds, less than call2
        final List<SVCallRecord> testCluster = Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff, SVTestUtils.call2);
        final SVCallRecord flattened = engine.getCollapser().apply(testCluster);
        Assert.assertEquals(flattened.getPositionA(), SVTestUtils.depthAndStuff.getPositionA());
        Assert.assertEquals(flattened.getPositionB(), SVTestUtils.depthAndStuff.getPositionB());
        //should have all the algs
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthOnly.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.call2.getAlgorithms()));
        //should have all the genotypes
        SVTestUtils.assertContainsAll(flattened.getGenotypes(), SVTestUtils.depthAndStuff.getGenotypes());
        SVTestUtils.assertContainsAll(flattened.getGenotypes(), SVTestUtils.depthOnly.getGenotypes());
        SVTestUtils.assertContainsAll(flattened.getGenotypes(), SVTestUtils.call2.getGenotypes());
        //TODO: add test for insertion cluster
    }

    @Test
    public void testClusterTogether() {
        Assert.assertTrue(engine.clusterTogether(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));
        Assert.assertFalse(engine.clusterTogether(SVTestUtils.depthOnly, SVTestUtils.inversion));
        Assert.assertFalse(engine.clusterTogether(SVTestUtils.call1, SVTestUtils.call2));
        Assert.assertTrue(engine.clusterTogether(SVTestUtils.call1, SVTestUtils.overlapsCall1));
    }

    @Test
    public void testGetClusteringIntervalEdge() {
        //edge case - end of contig
        Assert.assertTrue(engine.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) <= SVTestUtils.chr1Length);
    }

    @DataProvider(name = "recordPairs")
    public Object[][] recordPairs() {
        return new Object[][]{
                {SVTestUtils.call1, SVTestUtils.call1b, "depth-depth"},
                {SVTestUtils.call1, SVTestUtils.depthAndStuff2, "mixed"},
                {SVTestUtils.depthAndStuff, SVTestUtils.depthAndStuff2, "mixed"}
        };
    }

    @Test(dataProvider= "recordPairs")
    public void testGetMaxClusterableStartingPositionPairs(final SVCallRecord call1, final SVCallRecord call2, final String name) {
        //before we check calculations, make sure these actually cluster
        Assert.assertTrue(engine.clusterTogether(call1, call2), name);

        //max start for combined interval should be greater than the leftmost bound, and less than the nearest event end
        final int call1Position = engine.getMaxClusterableStartingPosition(call1);
        final int call2Position = engine.getMaxClusterableStartingPosition(call2);
        final int maxPosBoth = Math.max(call1Position, call2Position);
        Assert.assertTrue(maxPosBoth > call1.getPositionA(), name);
        Assert.assertTrue(maxPosBoth < call1.getPositionB(), name);
        Assert.assertTrue(maxPosBoth > call2.getPositionA(), name);
        Assert.assertTrue(maxPosBoth < call2.getPositionB(), name);

        //quantitative checks
        final int exactMaxPos1 = getMaxExactClusterablePosition(call1);
        final int exactMaxPos2 = getMaxExactClusterablePosition(call2);
        Assert.assertTrue(maxPosBoth >= exactMaxPos1, name);
        Assert.assertTrue(maxPosBoth >= exactMaxPos2, name);
        Assert.assertTrue(maxPosBoth >= Math.max(exactMaxPos1, exactMaxPos2), name);
    }

    final int getMaxExactClusterablePosition(final SVCallRecord record) {
        int pos = record.getPositionA() + 1;
        while (pos <= record.getPositionB()) {
            final SVCallRecord testRecord = new SVCallRecord("", record.getContigA(), pos, record.getStrandA(),
                    record.getContigB(), pos + record.getLength(), record.getStrandB(), record.getType(),
                    record.getLength(), record.getAlgorithms(), record.getAlleles(), record.getGenotypes());
            if (engine.clusterTogether(record, testRecord)) {
                pos++;
            } else {
                return pos - 1;
            }
        }
        throw new TestException("Max clusterable position not found");
    }

    @Test
    public void testGetClusteringIntervalLists() {
        //test overloaded function with List
        final List<SVCallRecord> pesrClusterList = new ArrayList<>();
        pesrClusterList.add(SVTestUtils.depthAndStuff);
        final int pesrCluster1 = engine.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + engine.getEvidenceParams().getWindow());
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + engine.getMixedParams().getWindow());
        //add an upstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff2);
        final int pesrCluster2 = engine.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster2 >= pesrCluster1);
        //add a downstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff3);
        final int pesrCluster3 = engine.getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + engine.getEvidenceParams().getWindow());
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + engine.getMixedParams().getWindow());
    }

    @Test
    public void testIsDepthOnlyCall() {
        Assert.assertTrue(SVTestUtils.call1.isDepthOnly());
        Assert.assertFalse(SVTestUtils.depthAndStuff.isDepthOnly());
        Assert.assertFalse(SVTestUtils.inversion.isDepthOnly());
    }

    @DataProvider(name = "clusterTogetherVaryPositionsProvider")
    public Object[][] clusterTogetherVaryPositionsProvider() {
        return new Object[][]{
                {500, 1001, 1001, 1502, false},  // abutting
                {500, 1001, 500, 1001, true},  // exactly equal
                {500, 1001, 600, 1101, true},  // barely meet reciprocal overlap
                {500, 1001, 601, 1102, false}, // call2 shifted slightly up
                {500, 1000, 600, 1101, false}, // call1 slightly larger
                {500, 501, 500, 501, true}, // tiny but equal
                {500, 501, 500, 502, false}, // tiny but call2 twice as big
                {500, 500, 500, 500, true}, // 0-length and equal
                {500, 500, 501, 501, false}, // 0-length and not equal
                {1, SVTestUtils.chr1Length, 1, SVTestUtils.chr1Length, true}, // really big
                {1, 10001, 1, 10001, true}, // left contig edge
                {SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, true}, // right contig edge
                {100000, 200001, 102000, 202001, true} // window test fail
        };
    }

    @Test(dataProvider= "clusterTogetherVaryPositionsProvider")
    public void testClusterTogetherVaryPositions(final int start1, final int end1, final int start2, final int end2, final boolean result) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start1, true,
                "chr1", end1, false,
                StructuralVariantType.DEL, end1 - start1 + 1, Lists.newArrayList("pesr"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap());
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", start2, true,
                "chr1", end2, false,
                StructuralVariantType.DEL, end2 - start2 + 1, Lists.newArrayList("depth"),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap());
        Assert.assertEquals(engine.clusterTogether(call1, call2), result);
    }

    @Test
    public void testClusterTogetherVaryTypes() {
        for (final StructuralVariantType type1 : StructuralVariantType.values()) {
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                    "chr1", 2001, false, type1,
                    1000, Lists.newArrayList("depth"),
                    Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
            for (final StructuralVariantType type2 : StructuralVariantType.values()) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, true,
                        "chr1", 2001, false, type2,
                        1000, Lists.newArrayList("depth"),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                // Should only cluster together if same type
                Assert.assertEquals(engine.clusterTogether(call1, call2), type1 == type2);
            }
        }
    }

    @Test
    public void testClusterTogetherVaryStrands() {
        final List<Boolean> bools = Lists.newArrayList(Boolean.TRUE, Boolean.FALSE);
        for (final Boolean strand1A : bools) {
            for (final Boolean strand1B : bools) {
                final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, strand1A,
                        "chr1", 2001, strand1B, StructuralVariantType.BND,
                        1000, Lists.newArrayList("depth"),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                for (final Boolean strand2A : bools) {
                    for (final Boolean strand2B : bools) {
                        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, strand2A,
                                "chr1", 2001, strand2B, StructuralVariantType.BND,
                                1000, Lists.newArrayList("depth"),
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                        // Should only cluster if strands match
                        Assert.assertEquals(engine.clusterTogether(call1, call2), strand1A == strand2A && strand1B == strand2B);
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryContigs() {
        final List<String> contigs = Lists.newArrayList("chr1", "chrX");
        for (final String contig1A : contigs) {
            for (final String contig1B : contigs) {
                final SVCallRecord call1 = new SVCallRecord("call1", contig1A, 1000, true,
                        contig1B, 2001, false, StructuralVariantType.BND,
                        1000, Lists.newArrayList("depth"),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                for (final String contig2A : contigs) {
                    for (final String contig2B : contigs) {
                        final SVCallRecord call2 = new SVCallRecord("call2", contig2A, 1000, true,
                                contig2B, 2001, false, StructuralVariantType.BND,
                                1000, Lists.newArrayList("depth"),
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                        // Should only cluster if contigs match
                        Assert.assertEquals(engine.clusterTogether(call1, call2), contig1A.equals(contig2A) && contig1B.equals(contig2B));
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryAlgorithms() {
        final List<List<String>> algorithmsList = Lists.newArrayList(
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, "PESR"),
                Arrays.asList("PESR")
        );
        for (final List<String> algorithms1 : algorithmsList) {
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                    "chr1", 2001, false, StructuralVariantType.DEL,
                    1000, algorithms1, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
            for (final List<String> algorithms2 : algorithmsList) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, true,
                        "chr1", 2001, false, StructuralVariantType.DEL,
                        1000, algorithms2, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
                // All combinations should cluster
                Assert.assertTrue(engine.clusterTogether(call1, call2));
            }
        }
    }

    @Test
    public void testClusterTogetherVaryParameters() {
        final SVClusterEngine<SVCallRecord> testEngine = SVTestUtils.getNewDefaultSingleLinkageEngine();
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                "chr1", 2001, false, StructuralVariantType.DEL,
                1000, Collections.singletonList("depth"), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1100, true,
                "chr1", 2101, false, StructuralVariantType.DEL,
                1000, Collections.singletonList("depth"), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
        // Cluster with default parameters
        Assert.assertTrue(testEngine.clusterTogether(call1, call2));
        final SVClusterEngine.ClusteringParameters exactMatchParameters = new SVClusterEngine.DepthClusteringParameters(1.0, 0, 1.0);
        testEngine.setDepthOnlyParams(exactMatchParameters);
        // Do not cluster requiring exact overlap
        Assert.assertFalse(testEngine.clusterTogether(call1, call2));
    }


    @DataProvider(name = "testAddVaryPositionsProvider")
    public Object[][] testAddVaryPositionsProvider() {
        return new Object[][]{
                {10000, 20001, 12000, 22001, 14000, 24001, SINGLE_LINKAGE, 1},
                {10000, 20001, 12000, 22001, 14001, 24002, SINGLE_LINKAGE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, SINGLE_LINKAGE, 3},
                {10000, 20001, 11000, 21001, 12000, 22001, MAX_CLIQUE, 1},
                {10000, 20001, 12000, 22001, 14000, 24001, MAX_CLIQUE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, MAX_CLIQUE, 3}
        };
    }

    @Test(dataProvider= "testAddVaryPositionsProvider")
    public void testAddVaryPositions(final int positionA1, final int positionB1,
                                     final int positionA2, final int positionB2,
                                     final int positionA3, final int positionB3,
                                     final LocatableClusterEngine.CLUSTERING_TYPE type,
                                     final int result) {
        final SVClusterEngine<SVCallRecord> engine;
        if (type == SINGLE_LINKAGE) {
            engine = SVTestUtils.getNewDefaultSingleLinkageEngine();
        } else if (type == MAX_CLIQUE) {
            engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        } else {
            throw new TestException("Unimplemented clustering type " + type.name());
        }
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", positionA1, true,
                "chr1", positionB1, false, StructuralVariantType.DEL,
                positionB1 - positionA1 + 1, Collections.singletonList("depth"),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
        final SVCallRecord call2 = new SVCallRecord("call1", "chr1", positionA2, true,
                "chr1", positionB2, false, StructuralVariantType.DEL,
                positionB2 - positionA2 + 1, Collections.singletonList("depth"),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
        final SVCallRecord call3 = new SVCallRecord("call1", "chr1", positionA3, true,
                "chr1", positionB3, false, StructuralVariantType.DEL,
                positionB3 - positionA3 + 1, Collections.singletonList("depth"),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap());
        engine.add(call1);
        engine.add(call2);
        engine.add(call3);
        Assert.assertEquals(engine.getOutput().size(), result);
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine<SVCallRecord> temp1 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        Assert.assertTrue(temp1.isEmpty());
        temp1.add(SVTestUtils.call1);
        Assert.assertFalse(temp1.isEmpty());
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertTrue(temp1.isEmpty());
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembership(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembership(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine<SVCallRecord> temp2 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.overlapsCall1);
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        //median of two items ends up being the second item here
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.overlapsCall1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.overlapsCall1.getPositionB());
        SVTestUtils.assertEqualsExceptMembership(output2.get(1), SVTestUtils.call4_chr10);

        //checking insensitivity to sample set overlap
        final SVClusterEngine<SVCallRecord> temp3 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecord> output3 = temp3.getOutput();
        Assert.assertEquals(output3.size(), 1);
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.call1.getPositionB());
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.sameBoundsSampleMismatch.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.sameBoundsSampleMismatch.getPositionB());
    }
}