package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import org.apache.commons.math3.util.FastMath;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.apache.commons.math3.util.FastMath.max;

class FilterLoss implements Comparable<FilterLoss> {
    final double inheritanceLoss;
    final double truthLoss;
    final String label;
    static final double truthWeight = MinGqVariantFilterBase.truthWeight;
    static final double inheritanceWeight = 1.0 - truthWeight;
    static final double targetPrecision = MinGqVariantFilterBase.targetPrecision;

    FilterLoss(final double inheritanceLoss, final double truthLoss, final String label) {
        this.inheritanceLoss = inheritanceLoss;
        this.truthLoss = truthLoss;
        this.label = label;
    }

    FilterLoss(final double inheritanceLoss, final double truthLoss) {
        this.inheritanceLoss = inheritanceLoss;
        this.truthLoss = truthLoss;
        this.label = null;
    }

    FilterLoss(final FilterSummary filterSummary) {
        this(1.0 - getInheritanceF1(filterSummary),1.0 - getTruthF1(filterSummary),
             filterSummary.label);
    }

    FilterLoss(final FilterLoss copyLoss) {
        this(copyLoss.inheritanceLoss, copyLoss.truthLoss, copyLoss.label);
    }

    FilterLoss(final FilterLoss copyLoss, final String label) {
        this(copyLoss.inheritanceLoss, copyLoss.truthLoss, label);
    }

    FilterLoss(final Collection<FilterSummary> filterSummaries) {
        this(
                filterSummaries.stream()
                        .map(FilterLoss::new)
                        .reduce(FilterLoss::add).orElse(EMPTY)
                        .divide(max((int) filterSummaries.stream()
                                        .filter(FilterSummary::isNotEmpty).count(),
                                1)),
                filterSummaries.stream()
                        .filter(FilterSummary::isNotEmpty)
                        .map(FilterLoss::new)
                        .map(FilterLoss::toString)
                        .collect(Collectors.joining("\n"))
                        + "\noverall: "
        );
    }

    protected static double getInheritanceF1(final FilterSummary filterSummary) {
        // calculate f1 score:
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        //     recall = numMendelian / numDiscoverableMendelian
        //     precision = numMendelian / numPassed
        //     -> f1 = 2.0 * numMendelian / (numNonRef + numPassed)
        // No data to filter? -> no loss -> f1 = 1.0
        if(filterSummary.numDiscoverable == 0) {
            // No way to score this if no Mendelian trios are discoverable
            return Double.NaN;
        } else if(filterSummary.numMendelian == 0) {
            return 0.0;  // Discovered nothing -> get bad score even if you passed nothing
        }
        final double precision = filterSummary.numMendelian / (double)filterSummary.numPassed;
        if(precision > targetPrecision) {
            final double recall = filterSummary.numMendelian / (double) filterSummary.numDiscoverable;
            final double f1 = 2.0 / (1.0 / recall + 1.0 / precision);
            return f1 * (1.0 - targetPrecision) + targetPrecision;
        } else {
            return precision; // just maximize precision
        }
    }

    protected static double getTruthF1(final FilterSummary filterSummary) {
        // NOTE: calls are compared against list of true variants and false variants
        // calculate f1 score:
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        //     recall = numTruePositive / (numTruePositive + numFalseNegative)
        //     precision = numTruePositive / (numTruePositive + numFalsePositive)
        //     -> f1 = 2.0 * numTruePositive / (2 * numTruePositive + numFalseNegative + numFalsePositive)
        // No data to filter? -> no loss -> f1 = 1.0
        final long numTrue = filterSummary.numTruePositive + filterSummary.numFalseNegative;
        final long numFalse = filterSummary.numTrueNegative + filterSummary.numFalsePositive;
        if(numFalse == 0) {
            if(numTrue == 0) {
                // no truth data whatsoever
                return Double.NaN;
            }
            // no "false" variants in truth set, have recall but no precision
            final double recall = filterSummary.numTruePositive / (double) numTrue;
            // if recall is above targetPrecision (i.e. presumably very good) try to get a good ratio of truePositives
            // to passed (i.e. pass the fewest number of unknown variants)
            final double trueToPassed = filterSummary.numTruePositive / (double) filterSummary.numPassed;
            return recall <= targetPrecision ? recall : targetPrecision + (1.0 - targetPrecision) * trueToPassed;
        } else if(numTrue == 0) {
            // no "true" variants in truth set, but can optimize for not including garbage
            final double negativeRecall = filterSummary.numTrueNegative / (double) numFalse;
            // if "negative recall" is good enough, try to pass as much as possible
            final double passRatio = filterSummary.numPassed / (double)(filterSummary.numPassed + numFalse);
            return negativeRecall <= targetPrecision ? negativeRecall : targetPrecision + (1.0 - targetPrecision) * passRatio;
        }
        final long numClassifiedPositive = filterSummary.numTruePositive + filterSummary.numFalsePositive;
        if(numClassifiedPositive == 0) {
            // we know there is truth data, so this is just a crappy classifier
            return 0.0;
        }
        // this is guaranteed to be finite now
        final double precision = filterSummary.numTruePositive / (double) numClassifiedPositive;
        if(precision <= targetPrecision) {
            return precision;  // haven't met target precision, just optimize that
        }
        final double recall = filterSummary.numTruePositive / (double)numTrue;
        final double f1 = 2.0 / (1.0 / recall + 1.0 / precision);
        return f1 * (1.0 - targetPrecision) + targetPrecision;
    }

    FilterLoss(final float[] pSampleVariantIsGood, final boolean[] sampleVariantIsGood) {
        final double balancedLoss = IntStream.range(0, pSampleVariantIsGood.length)
            .mapToDouble(
                idx -> -FastMath.log(
                    sampleVariantIsGood[idx] ?
                        (double) pSampleVariantIsGood[idx] :
                        1.0 - (double) pSampleVariantIsGood[idx]
                )
            )
            .sum();
        this.inheritanceLoss = balancedLoss;
        this.truthLoss = balancedLoss;
        this.label = null;
    }

    static FilterLoss add(final FilterLoss lossA, final FilterLoss lossB) {
        return new FilterLoss(addLosses(lossA.inheritanceLoss, lossB.inheritanceLoss),
                              addLosses(lossA.truthLoss, lossB.truthLoss),
                              lossA.label == null ? lossB.label : lossA.label);
    }

    static private double addLosses(final double lossA, final double lossB) {
        if(Double.isFinite(lossA)) {
            if(Double.isFinite(lossB)) {
                return lossA + lossB;
            } else {
                return lossA;
            }
        } else {
            return lossB;
        }
    }

    FilterLoss divide(final double scale) {
        return new FilterLoss(inheritanceLoss / scale, truthLoss / scale, label);
    }

    public boolean ge(final FilterLoss other) {
        return toDouble() >= other.toDouble();
    }

    public boolean gt(final FilterLoss other) {
        return toDouble() > other.toDouble();
    }

    public boolean le(final FilterLoss other) {
        return toDouble() <= other.toDouble();
    }

    public boolean lt(final FilterLoss other) {
        return toDouble() < other.toDouble();
    }

    @Override
    public int compareTo(final @NotNull FilterLoss other) {
        if (this.lt(other)) {
            return -1;
        } else if (this.gt(other)) {
            return 1;
        } else {
            return 0;
        }
    }

    private static String getDescription(final double inheritanceLoss, final double truthLoss) {
        return "{inherit:" + inheritanceLoss + ",truth:" + truthLoss + "}";
    }

    @Override
    public String toString() {
        return label == null ?
                getDescription(inheritanceLoss, truthLoss) :
                label + ": " + getDescription(inheritanceLoss, truthLoss);
    }

    double toDouble() {
        final double weightedInheritanceLoss = Double.isFinite(inheritanceLoss) ?
            inheritanceWeight * inheritanceLoss :
            0.0;
        final double weightedTruthLoss = Double.isFinite(truthLoss) ?
            truthWeight * truthLoss :
            0.0;
        return weightedInheritanceLoss + weightedTruthLoss;
    }

    float toFloat() {
        return (float) toDouble();
    }

    static final FilterLoss EMPTY = new FilterLoss(FilterSummary.EMPTY);
}
