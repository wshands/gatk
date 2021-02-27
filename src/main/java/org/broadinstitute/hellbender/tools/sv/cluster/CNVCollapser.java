package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Set;
import java.util.stream.Collectors;

public class CNVCollapser extends SVCollapser<SVCallRecord> {

    public CNVCollapser(final BreakpointSummaryStrategy strategy) {
        super(strategy);
    }

    @Override
    public SVCallRecord collapse(Collection<SVCallRecord> items) {
        return collapseToBaseRecord(items);
    }

    @Override
    protected Genotype collapseSampleGenotypes(final Collection<Genotype> genotypes) {
        Utils.nonEmpty(genotypes);
        final String sampleName = genotypes.iterator().next().getSampleName();
        if (!genotypes.stream().allMatch(g -> g.getSampleName().equals(sampleName))) {
            throw new IllegalArgumentException("This method expects a list of genotypes from the same sample, " +
                    "but not all input genotypes represent sample " + sampleName + ".");
        }
        final Set<Integer> copyNumbers = genotypes.stream()
                .filter(g -> g.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT))
                .map(g -> Integer.parseInt(g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT).toString()))
                .collect(Collectors.toSet());
        final GenotypeBuilder gb = new GenotypeBuilder(genotypes.iterator().next());
        if (copyNumbers.isEmpty()) {
            return gb.make();
        } else if (copyNumbers.size() == 1) {
            //For now just make sure genotypes have the same copy number -- qualities will be recalculated elsewhere
            gb.noAttributes();
            gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumbers.iterator().next());
            return gb.make();
        } else {
            throw new IllegalArgumentException("This method will only merge genotypes with the same copy number. Found " + copyNumbers.size() + " different copy numbers.");
        }
    }
}
