package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class SVClusterEngineArgumentsCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    protected static final String BASE_OVERLAP_FRACTION_NAME = "-overlap-fraction";
    protected static final String BASE_BREAKEND_WINDOW_NAME = "-breakend-window";

    protected static final String BASE_OVERLAP_FRACTION_DOC = "variant reciprocal overlap fraction";
    protected static final String BASE_BREAKEND_WINDOW_DOC = " window size for breakend proximity";

    private static final String DEPTH_OVERLAP_FRACTION_NAME = "depth" + BASE_OVERLAP_FRACTION_NAME;
    private static final String MIXED_OVERLAP_FRACTION_NAME = "mixed" + BASE_OVERLAP_FRACTION_NAME;
    private static final String PESR_OVERLAP_FRACTION_NAME = "pesr" + BASE_OVERLAP_FRACTION_NAME;

    private static final String DEPTH_BREAKEND_WINDOW_NAME = "depth" + BASE_BREAKEND_WINDOW_NAME;
    private static final String MIXED_BREAKEND_WINDOW_NAME = "mixed" + BASE_BREAKEND_WINDOW_NAME;
    private static final String PESR_BREAKEND_WINDOW_NAME = "pesr" + BASE_BREAKEND_WINDOW_NAME;

    @Argument(fullName = DEPTH_OVERLAP_FRACTION_NAME,
            doc="Depth/Depth" + BASE_OVERLAP_FRACTION_DOC, optional=true)
    public double depthOverlapFraction = 0.8;

    @Argument(fullName = MIXED_OVERLAP_FRACTION_NAME,
            doc="PESR/Depth" + BASE_OVERLAP_FRACTION_DOC, optional=true)
    public double mixedOverlapFraction = 0.5;

    @Argument(fullName = PESR_OVERLAP_FRACTION_NAME,
            doc="PESR/PESR" + BASE_OVERLAP_FRACTION_DOC, optional=true)
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

    public final SVClusterEngine.DepthClusteringParameters getDepthParameters() {
        return new SVClusterEngine.DepthClusteringParameters(depthOverlapFraction, depthBreakendWindow);
    }

    public final SVClusterEngine.MixedClusteringParameters getMixedParameters() {
        return new SVClusterEngine.MixedClusteringParameters(mixedOverlapFraction, mixedBreakendWindow);
    }

    public final SVClusterEngine.EvidenceClusteringParameters getPESRParameters() {
        return new SVClusterEngine.EvidenceClusteringParameters(pesrOverlapFraction, pesrBreakendWindow);
    }
}
