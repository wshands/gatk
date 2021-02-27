package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.Collection;

public class SVPreprocessingRecordCollapser extends SVPreprocessingCollapser<SVCallRecord> {

    public SVPreprocessingRecordCollapser(final BreakpointSummaryStrategy strategy) {
        super(strategy);
    }

    @Override
    public SVCallRecord collapse(final Collection<SVCallRecord> items) {
        return collapseToBaseRecord(items);
    }
}
