package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class SVClusterEngineTest {

    private final SVClusterEngine<SVCallRecord> engine = SVTestUtils.defaultEngine;

    @BeforeTest
    public void initializeClusterEngine() {
        engine.add(SVTestUtils.call1);
    }

    @Test
    public void testFlattenCluster() {
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
    public void testGetClusteringInterval() {
        Assert.assertTrue(engine.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) <= SVTestUtils.chr1Length);

        final int littleClusterPosition = engine.getMaxClusterableStartingPosition(SVTestUtils.call1);
        //before we check calculations, make sure these actually cluster
        Assert.assertTrue(engine.clusterTogether(SVTestUtils.call1, SVTestUtils.call1b));
        final int callPosition = engine.getMaxClusterableStartingPosition(SVTestUtils.call1b);
        final int totalPosition = Math.max(callPosition, littleClusterPosition);
        //max start for combined interval should be greater than the leftmost bound, and less than the nearest event end
        Assert.assertTrue(totalPosition > SVTestUtils.call1.getPositionA());
        Assert.assertTrue(totalPosition < SVTestUtils.call1.getPositionB());
        //quantitative checks
        final double depthReciprocalOverlap = Math.max(engine.getDepthOnlyParams().getReciprocalOverlap(), engine.getMixedParams().getReciprocalOverlap());
        final int depthPadding = Math.max(engine.getDepthOnlyParams().getPadding(), engine.getMixedParams().getPadding());
        final int depthWindow = Math.max(engine.getDepthOnlyParams().getWindow(), engine.getMixedParams().getWindow());
        final SimpleInterval paddedInterval = new SimpleInterval(SVTestUtils.call1b.getContigA(), SVTestUtils.call1b.getPositionA(), SVTestUtils.call1b.getPositionB()).expandWithinContig(depthPadding, SVTestUtils.dict);
        Assert.assertTrue(totalPosition >= paddedInterval.getStart() + (1.0 - depthReciprocalOverlap) * paddedInterval.getLengthOnReference());
        Assert.assertTrue(totalPosition >= SVTestUtils.call1b.getPositionA() + depthWindow);

        //checks for breakend calls, which have different bounds
        final int pesrWindow = Math.max(engine.getEvidenceParams().getWindow(), engine.getMixedParams().getWindow());
        final int pesrCluster1 = engine.getMaxClusterableStartingPosition(SVTestUtils.depthAndStuff);
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + pesrWindow);
        //add an upstream variant
        final int pesrCluster2 = Math.max(engine.getMaxClusterableStartingPosition(SVTestUtils.depthAndStuff2), pesrCluster1);
        Assert.assertTrue(pesrCluster2 >= pesrCluster1);
        //add a downstream variant
        final int pesrCluster3 = Math.max(engine.getMaxClusterableStartingPosition(SVTestUtils.depthAndStuff3), pesrCluster1);
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + pesrWindow);
    }

    @Test
    public void testIsDepthOnlyCall() {
        Assert.assertTrue(SVClusterEngine.isDepthOnlyCall(SVTestUtils.depthOnly));
        Assert.assertFalse(SVClusterEngine.isDepthOnlyCall(SVTestUtils.depthAndStuff));
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine temp1 = SVTestUtils.getNewDefaultEngine();
        temp1.add(SVTestUtils.call1);
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEquals(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEquals(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine temp2 = SVTestUtils.getNewDefaultEngine();
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.overlapsCall1);
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        //median of two items ends up being the second item here
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.overlapsCall1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.overlapsCall1.getPositionB());
        SVTestUtils.assertEquals(output2.get(1), SVTestUtils.call4_chr10);

        //checking insensitivity to sample set overlap
        final SVClusterEngine temp3 = SVTestUtils.getNewDefaultEngine();
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