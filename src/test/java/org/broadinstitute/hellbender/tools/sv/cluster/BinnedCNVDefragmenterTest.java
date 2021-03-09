package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BinnedCNVDefragmenterTest {

    private static final double paddingFraction = 0.25;
    private static final double sampleOverlap = 0.9;
    private static final CNVDefragmenter defaultDefragmenter = new CNVDefragmenter(SVTestUtils.dict, paddingFraction, sampleOverlap);
    private static final CNVDefragmenter singleSampleDefragmenter = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0, SVTestUtils.targetIntervals);

    @Test
    public void testFlattenCluster() {
        final SVCallRecord call1FlattenedDefault = defaultDefragmenter.getCollapser().apply(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(SVTestUtils.call1, call1FlattenedDefault);

        final SVCallRecord call1FlattenedSingleSample = singleSampleDefragmenter.getCollapser().apply(Collections.singletonList(SVTestUtils.call1));
        Assert.assertEquals(call1FlattenedSingleSample, call1FlattenedDefault);

        final SVCallRecord sameBoundsThreeSamples = singleSampleDefragmenter.getCollapser().apply(Arrays.asList(SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch));
        Assert.assertEquals(sameBoundsThreeSamples.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(sameBoundsThreeSamples.getPositionB(), SVTestUtils.call1.getPositionB());
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(0).sameGenotype(SVTestUtils.sample1));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(1).sameGenotype(SVTestUtils.sample2));
        Assert.assertTrue(sameBoundsThreeSamples.getGenotypes().get(2).sameGenotype(SVTestUtils.sample3));

        final SVCallRecord overlapping = singleSampleDefragmenter.getCollapser().apply(Arrays.asList(SVTestUtils.call1, SVTestUtils.call2));
        Assert.assertEquals(overlapping.getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(overlapping.getPositionB(), SVTestUtils.call2.getPositionB());
    }

    @DataProvider
    public Object[][] clusterTogetherInputsDefault() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true},
                {SVTestUtils.call1, SVTestUtils.call2, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, false},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false}
        };
    }

    @DataProvider
    public Object[][] clusterTogetherInputsSingleSample() {
        return new Object[][] {
                {SVTestUtils.call1, SVTestUtils.call1, true},
                {SVTestUtils.call1, SVTestUtils.call2, true},  //overlapping, same samples
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1, SVTestUtils.sameBoundsSampleMismatch, true},
                {SVTestUtils.call1, SVTestUtils.nonDepthOnly, false},
                {SVTestUtils.call1_CN1, SVTestUtils.call2_CN0, false}  //overlapping, but different copy number
        };
    }

    @Test(dataProvider = "clusterTogetherInputsDefault")
    public void testClusterTogetherDefault(final SVCallRecord call1, final SVCallRecord call2, final boolean expectedResult) {
        Assert.assertEquals(defaultDefragmenter.clusterTogether(call1, call2), expectedResult);
    }

    @Test(dataProvider = "clusterTogetherInputsSingleSample")
    public void testClusterTogetherSingleSample(final SVCallRecord call1, final SVCallRecord call2, final boolean expectedResult) {
        Assert.assertEquals(singleSampleDefragmenter.clusterTogether(call1, call2), expectedResult);
    }

    @Test
    public void testGetClusteringInterval() {
        Assert.assertEquals(defaultDefragmenter.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall), SVTestUtils.chr1Length);
        Assert.assertTrue(singleSampleDefragmenter.getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) == SVTestUtils.chr1Length);  //will be less than chr1length if target intervals are smaller than chr1

        final int totalPosition = defaultDefragmenter.getMaxClusterableStartingPosition(SVTestUtils.call2);
        //padding is added to the input call
        //Assert.assertEquals(totalPosition, SVTestUtils.call2.getPositionA()+(int)Math.round(SVTestUtils.length*defaultDefragmenter.getPaddingFraction()));
        // TODO this needs to be reworked
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final CNVDefragmenter temp1 = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        temp1.add(SVTestUtils.call1);
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.getOutput(); //flushes all clusters
        Assert.assertEquals(output1.size(), 2);
        Assert.assertEquals(SVTestUtils.call1, output1.get(0));
        Assert.assertEquals(SVTestUtils.call3, output1.get(1));

        final CNVDefragmenter temp2 = new BinnedCNVDefragmenter(SVTestUtils.dict, paddingFraction, 0.8, SVTestUtils.targetIntervals);
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.call2);  //should overlap after padding
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.getOutput();
        Assert.assertEquals(output2.size(), 2);
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call2.getPositionB());
        Assert.assertEquals(output2.get(1), SVTestUtils.call4_chr10);

        //cohort case, checking sample set overlap
        final CNVDefragmenter temp3 = new CNVDefragmenter(SVTestUtils.dict);
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecord> output3 = temp3.getOutput();
        Assert.assertEquals(output3.size(), 2);
    }
}