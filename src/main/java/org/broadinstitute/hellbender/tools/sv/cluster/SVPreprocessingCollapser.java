package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.*;

public abstract class SVPreprocessingCollapser<T extends SVCallRecord> extends SVCollapser<T> {

    public SVPreprocessingCollapser(final BreakpointSummaryStrategy strategy) {
        super(strategy);
    }

    @Override
    protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes) {
        final List<Allele> alleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
        final GenotypeBuilder builder = new GenotypeBuilder(genotypes.iterator().next());
        builder.alleles(alleles);
        builder.noAttributes();
        builder.attributes(collapseAttributes(genotypes));
        return builder.make();
    }

    protected Map<String, Object> collapseAttributes(final Collection<Genotype> genotypes) {
        if (genotypes.stream().anyMatch(g -> SVCallRecord.isRawCall(g) || SVCallRecord.isCarrier(g))) {
            return Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_TRUE);
        } else {
            return Collections.singletonMap(GATKSVVCFConstants.RAW_CALL_ATTRIBUTE, GATKSVVCFConstants.RAW_CALL_ATTRIBUTE_FALSE);
        }
    }

}
