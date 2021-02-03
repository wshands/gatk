package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.sv.PairedEndAndSplitReadEvidenceCollection.LocusDepth;

@BetaFeature
@CommandLineProgramProperties(
        summary = "Prints allele counts.",
        oneLineSummary = "Prints allele counts.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
public class PrintAlleleCounts extends FeatureWalker<LocusDepth> {
    public static final String ALLELE_COUNT_INPUT_ARGUMENT_SHORT_NAME = "F";
    public static final String ALLELE_COUNT_INPUT_ARGUMENT_LONG_NAME = "allele-count-vcf";

    @Argument(shortName = ALLELE_COUNT_INPUT_ARGUMENT_SHORT_NAME,
            fullName = ALLELE_COUNT_INPUT_ARGUMENT_LONG_NAME,
            doc = "Input VCF of SNPs marking loci for allele counts")
    public String alleleCountInputFilename;

    @Override
    protected boolean isAcceptableFeatureType( final Class<? extends Feature> featureType ) {
        return featureType == LocusDepth.class;
    }

    @Override
    public void apply( final LocusDepth feature,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext, FeatureContext featureContext ) {
        System.out.println(feature.toString());
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return new GATKPath(alleleCountInputFilename);
    }
}
