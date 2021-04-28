package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

class FilterSummary {
    final int minGq;
    final long numDiscoverable;
    final long numPassed;
    final long numMendelian;
    final long numTruePositive;
    final long numFalsePositive;
    final long numFalseNegative;
    final long numTrueNegative;
    final String label;

    FilterSummary(final int minGq,
                  final long numDiscoverable, final long numPassed, final long numMendelian,
                  final long numTruePositive, final long numFalsePositive, final long numFalseNegative,
                  final long numTrueNegative, final String label) {
        this.minGq = minGq;
        this.numDiscoverable = numDiscoverable;
        this.numPassed = numPassed;
        this.numMendelian = numMendelian;
        this.numTruePositive = numTruePositive;
        this.numFalsePositive = numFalsePositive;
        this.numFalseNegative = numFalseNegative;
        this.numTrueNegative = numTrueNegative;
        this.label = label;
    }

    FilterSummary(final FilterSummary other) {
        this.minGq = other.minGq;
        this.numDiscoverable = other.numDiscoverable;
        this.numPassed = other.numPassed;
        this.numMendelian = other.numMendelian;
        this.numTruePositive = other.numTruePositive;
        this.numFalsePositive = other.numFalsePositive;
        this.numFalseNegative = other.numFalseNegative;
        this.numTrueNegative = other.numTrueNegative;
        this.label = other.label;
    }

    static final FilterSummary EMPTY = new FilterSummary(
            Integer.MIN_VALUE, 0L, 0L, 0L, 0L,
            0L, 0L, 0L, null
    );

    boolean hasInheritanceData() { return this.numDiscoverable > 0; }

    boolean hasOverlapData() {
        return this.numFalsePositive + this.numTruePositive + this.numFalseNegative + this.numTrueNegative > 0;
    }

    boolean isEmpty() {
        return !isNotEmpty();
    }

    boolean isNotEmpty() {
        return hasInheritanceData() || hasOverlapData();
    }

    FilterSummary setMinGq(final int newMinGq) {
        return new FilterSummary(newMinGq, numDiscoverable, numPassed, numMendelian,
                numTruePositive, numFalsePositive, numFalseNegative, numTrueNegative, label);
    }

    FilterSummary shiftMinGq(final int minGqShift) {
        return setMinGq(minGq + minGqShift);
    }

    FilterSummary add(final FilterSummary other) {
        return new FilterSummary(
                minGq,
                numDiscoverable + other.numDiscoverable, numPassed + other.numPassed,
                numMendelian + other.numMendelian, numTruePositive + other.numTruePositive,
                numFalsePositive + other.numFalsePositive, numFalseNegative + other.numFalseNegative,
                numTrueNegative + other.numTrueNegative, label == null ? other.label : label
        );
    }

    FilterSummary subtract(final FilterSummary other) {
        return new FilterSummary(
                minGq,
                numDiscoverable - other.numDiscoverable, numPassed - other.numPassed,
                numMendelian - other.numMendelian, numTruePositive - other.numTruePositive,
                numFalsePositive - other.numFalsePositive, numFalseNegative - other.numFalseNegative,
                numTrueNegative - other.numTrueNegative, label == null ? other.label : label
        );
    }

    @Override
    public String toString() {
        return label + ":{minGq:" + minGq + ",  numDiscoverable:" + numDiscoverable + ", numPassed:" + numPassed
                + ", numMendelian:" + numMendelian + ", numTruePositive:" + numTruePositive
                + ", numFalsePositive:" + numFalsePositive + ", numFalseNegative: " + numFalseNegative
                + ", numTrueNegative:" + numTrueNegative + "}";
    }
}
