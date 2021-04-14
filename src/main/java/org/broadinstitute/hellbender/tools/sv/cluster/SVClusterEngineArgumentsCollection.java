package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class SVClusterEngineArgumentsCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_NAME = "-overlap-fraction";
    protected static final String BASE_BREAKEND_WINDOW_NAME = "-breakend-window";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_NAME = "-sample-overlap";

    protected static final String BASE_INTERVAL_OVERLAP_FRACTION_DOC = " interval reciprocal overlap fraction";
    protected static final String BASE_BREAKEND_WINDOW_DOC = " window size for breakend proximity";
    protected static final String BASE_SAMPLE_OVERLAP_FRACTION_DOC = " shared sample overlap fraction";

    public static final String DEPTH_INTERVAL_OVERLAP_FRACTION_NAME = "depth" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String MIXED_INTERVAL_OVERLAP_FRACTION_NAME = "mixed" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;
    public static final String PESR_INTERVAL_OVERLAP_FRACTION_NAME = "pesr" + BASE_INTERVAL_OVERLAP_FRACTION_NAME;

    public static final String DEPTH_BREAKEND_WINDOW_NAME = "depth" + BASE_BREAKEND_WINDOW_NAME;
    public static final String MIXED_BREAKEND_WINDOW_NAME = "mixed" + BASE_BREAKEND_WINDOW_NAME;
    public static final String PESR_BREAKEND_WINDOW_NAME = "pesr" + BASE_BREAKEND_WINDOW_NAME;

    public static final String DEPTH_SAMPLE_OVERLAP_FRACTION_NAME = "depth" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;
    public static final String MIXED_SAMPLE_OVERLAP_FRACTION_NAME = "mixed" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;
    public static final String PESR_SAMPLE_OVERLAP_FRACTION_NAME = "pesr" + BASE_SAMPLE_OVERLAP_FRACTION_NAME;

    @Argument(fullName = DEPTH_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true)
    public double depthOverlapFraction = 0.8;

    @Argument(fullName = MIXED_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/Depth" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true)
    public double mixedOverlapFraction = 0.5;

    @Argument(fullName = PESR_INTERVAL_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_INTERVAL_OVERLAP_FRACTION_DOC, optional=true)
    public double pesrOverlapFraction = 0.5;

    @Argument(fullName = DEPTH_BREAKEND_WINDOW_NAME,
            doc="Depth/Depth" + BASE_BREAKEND_WINDOW_DOC, optional=true)
    public int depthBreakendWindow = 0;

    @Argument(fullName = MIXED_BREAKEND_WINDOW_NAME,
            doc="Depth/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true)
    public int mixedBreakendWindow = 1000;

    @Argument(fullName = PESR_BREAKEND_WINDOW_NAME,
            doc="PESR/PESR" + BASE_BREAKEND_WINDOW_DOC, optional=true)
    public int pesrBreakendWindow = 500;

    @Argument(fullName = DEPTH_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true)
    public double depthSampleOverlapFraction = 0.8;

    @Argument(fullName = MIXED_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="Depth/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true)
    public double mixedSampleOverlapFraction = 0.8;

    @Argument(fullName = PESR_SAMPLE_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_SAMPLE_OVERLAP_FRACTION_DOC, optional=true)
    public double pesrSampleOverlapFraction = 0.8;

    public final SVClusterEngine.DepthClusteringParameters getDepthParameters() {
        return new SVClusterEngine.DepthClusteringParameters(depthOverlapFraction, depthBreakendWindow, depthSampleOverlapFraction);
    }

    public final SVClusterEngine.MixedClusteringParameters getMixedParameters() {
        return new SVClusterEngine.MixedClusteringParameters(mixedOverlapFraction, mixedBreakendWindow, mixedSampleOverlapFraction);
    }

    public final SVClusterEngine.EvidenceClusteringParameters getPESRParameters() {
        return new SVClusterEngine.EvidenceClusteringParameters(pesrOverlapFraction, pesrBreakendWindow, pesrSampleOverlapFraction);
    }
}
