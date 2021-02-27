package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecordWithEvidence;

import java.util.Collection;

public class SVPreprocessingRecordWithEvidenceCollapser extends SVPreprocessingCollapser<SVCallRecordWithEvidence> {

    public SVPreprocessingRecordWithEvidenceCollapser(final BreakpointSummaryStrategy strategy) {
        super(strategy);
    }

    @Override
    public SVCallRecordWithEvidence collapse(final Collection<SVCallRecordWithEvidence> items) {
        final SVCallRecordWithEvidence example = items.iterator().next();
        return new SVCallRecordWithEvidence(
                collapseToBaseRecord(items),
                example.getStartSplitReadSites(),
                example.getEndSplitReadSites(),
                example.getDiscordantPairs(),
                example.getCopyNumberDistribution()
        );
    }
}
