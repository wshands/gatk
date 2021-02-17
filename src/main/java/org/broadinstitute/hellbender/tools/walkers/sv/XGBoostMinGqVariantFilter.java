package org.broadinstitute.hellbender.tools.walkers.sv;

import ml.dmlc.xgboost4j.java.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

@CommandLineProgramProperties(
        summary = "Extract matrix of properties for each variant. Also extract, num_variants x num_trios x 3 tensors of" +
                "allele count and genotype quality. These data will be used to train a variant filter based on min GQ" +
                "(and stratified by other variant properties) that maximizes the admission of variants with Mendelian" +
                "inheritance pattern while omitting non-Mendelian variants." +
                "Derived class must implement abstract method train_filter()",
        oneLineSummary = "Extract data for training min GQ variant filter from .vcf and .ped file with trios.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class XGBoostMinGqVariantFilter extends MinGqVariantFilterBase {
    @Argument(fullName="max-training-rounds", shortName="tr", doc="Maximum number of rounds of training", optional=true, minValue=1)
    public int maxTrainingRounds = 100;
    @Argument(fullName="early-stopping-rounds", shortName="e", doc="Stop training if no improvement is made in validation set for this many rounds. Set <= 0 to disable.", optional=true)
    public int earlyStoppingRounds = 10;
    @Argument(fullName="initial-min-gq-quantile", shortName="q", doc="Initial guess for min GQ, as a quantile of gq in variants of trios.", optional=true)
    public double initialMinGqQuantile = 0.05;
    @Argument(fullName="learning-rate", shortName="lr", doc="Learning rate for xgboost", optional=true)
    public double eta = 1.0;
    @Argument(fullName="max-depth", shortName="d", doc="Max depth of boosted decision tree", optional=true, minValue=1)
    public int maxDepth = 6;
    @Argument(fullName="gamma", doc="Regularization factor for xgboost", optional=true)
    public double gamma = 1.0e-9;
    @Argument(fullName="subsample", doc="Proportion of data selected for each tree", optional=true)
    public double subsample = 0.9;
    @Argument(fullName="colsample-by-tree", doc="Proportion of columns selected for each tree", optional=true)
    public double colsampleByTree = 0.9;
    @Argument(fullName="colsample-by-level", doc="Proportion of columns selected for each level of each tree", optional=true)
    public double colsampleByLevel = 1.0;
    @Argument(fullName="min-child-weight", doc="Proportion of columns selected for each level of each tree", optional=true, minValue=0.0)
    public double minChildWeight = 1.0;
    @Argument(fullName="prediction-scale-factor", doc="Scale factor for raw predictions from xgboost", optional=true)
    public double predictionScaleFactor = 1.0;
    @Argument(fullName="max-minibatch-size", shortName="mbs", doc="Max number of sample variants to include in any mini batch",
              optional=true)
    public int maxMiniBatchSize = Integer.MAX_VALUE; // 100000;

    private Booster booster = null;
    private static final String TRAIN_MAT_KEY = "train";
    private static final String VALIDATION_MAT_KEY = "validation";


    protected float[] getRowMajorSampleProperties(int[] variantIndices) {
        // Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
        final int numRows = (int) getNumTrainableSampleVariants(variantIndices);

        final float[] rowMajorVariantProperties = new float[numRows * getNumProperties()];
        // Loop over variants and filterable samples. Store properties for each sample in a single flat array
        int flatIndex = 0;
        for(final int variantIndex : variantIndices) {
            for(int sampleIndex = 0; sampleIndex < getNumSamples(); ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                flatIndex = propertiesTable.copyPropertiesRow(
                    rowMajorVariantProperties, flatIndex, variantIndex, sampleIndex, needsNormalizedProperties()
                );
            }
        }

        return rowMajorVariantProperties;
    }

    private DMatrix getDMatrix(final int[] variantIndices) {
        final float[] propertiesArr = getRowMajorSampleProperties(variantIndices);
        final boolean[] sampleVariantTruth = getSampleVariantTruth(variantIndices);
        final int numRows = sampleVariantTruth.length;
        return getDMatrix(propertiesArr, sampleVariantTruth, numRows, getNumProperties());
    }

    private DMatrix getDMatrix(final float[] sampleVariantProperties) {
        final boolean[] sampleIsGood = new boolean[1];
        return getDMatrix(sampleVariantProperties, sampleIsGood, 1, sampleVariantProperties.length);
    }

    private DMatrix getDMatrix(final float[] propertiesArr, final boolean[] sampleVariantTruth,
                               final int numRows, final int numColumns) {
        try {
            for(final float val : propertiesArr) {
                if(!Float.isFinite(val)) {
                    throw new GATKException("rowMajorVariantProperties contains a non-finite value (" + val + ")");
                }
            }
            final DMatrix dMatrix = new DMatrix(
                    propertiesArr, numRows, numColumns, Float.NaN
            );
            // Set baseline (initial prediction for min GQ)
            final float[] baseline = new float[numRows];
            final float[] weights = new float[numRows];
            final float[] labels = new float[numRows];
            for(int idx = 0; idx < numRows; ++idx) {
                baseline[idx] = 0F;
                weights[idx] = 1F;
                labels[idx] = sampleVariantTruth[idx] ? 1F : 0F;
            }
            dMatrix.setLabel(labels);
            dMatrix.setBaseMargin(baseline);
            dMatrix.setWeight(weights);
            return dMatrix;
        }
        catch(XGBoostError xgBoostError) {
            throw new GATKException("Error forming DMatrix", xgBoostError);
        }
    }

    private Map<String, Object> getXgboostParams() {
        return new HashMap<String, Object>() {
            private static final long serialVersionUID = 0L;
            {
                put("eta", eta);
                put("max_depth", maxDepth);
                put("gamma", gamma);
                put("subsample", subsample);
                put("colsample_bytree", colsampleByTree);
                put("colsample_bylevel", colsampleByLevel);
                put("min_child_weight", minChildWeight);
                put("validate_parameters", true);
                put("objective", "binary:logistic");
                put("eval_metric", "logloss");
            }
        };
    }

    private float sigma(final float predict) {
        return 1.0F / (1.0F + (float)FastMath.exp(-predictionScaleFactor * predict));
    }

    private class DataSubset {
        final String name;
        final int[] variantIndices;
        final boolean isTrainingSet;
        final int maxMiniBatchSize;
        final int numFilterableSampleVariants;
        final int[] miniBatchSplits;

        final List<Float> scores;
        private int bestScoreInd;
        private float bestScore;

        final DMatrix dMatrix;
        final float[] pSampleVariantGood;
        final boolean[] sampleVariantTruth;
        final float[] d1Loss;
        final float[] d2Loss;
        final float[][] lossDerivatives;
        final double bestPossibleLoss;

        DataSubset(final String name, final int[] variantIndices, final boolean isTrainingSet, final int maxMiniBatchSize) {
            this.name = name;
            this.variantIndices = variantIndices;
            this.isTrainingSet = isTrainingSet;
            this.maxMiniBatchSize = maxMiniBatchSize;
            numFilterableSampleVariants = (int) getNumTrainableSampleVariants(variantIndices);
            miniBatchSplits = getMiniBatchSplits(variantIndices, maxMiniBatchSize, numFilterableSampleVariants);

            scores = new ArrayList<>();
            bestScore = Float.POSITIVE_INFINITY;
            bestScoreInd = 0;

            this.dMatrix = miniBatchSplits == null ? getDMatrix(variantIndices) : null;
            System.out.format("DataSubset %s processed in %d mini batches\n", name, getNumMiniBatches());

            // pre-allocate results if they are needed and we're not doing mini-batches
            final boolean preAllocateResults = miniBatchSplits == null;
            final boolean preAllocateDerivs = preAllocateResults && isTrainingSet;
            pSampleVariantGood = preAllocateResults ? new float[numFilterableSampleVariants] : null;
            d1Loss = preAllocateDerivs ? new float[numFilterableSampleVariants] : null;
            d2Loss = preAllocateDerivs ? new float[numFilterableSampleVariants] : null;
            lossDerivatives = new float[][] {d1Loss, d2Loss};

            sampleVariantTruth = getSampleVariantTruth(variantIndices);
            final float[] tempTruthProbs = new float[numFilterableSampleVariants];
            for(int idx = 0; idx < numFilterableSampleVariants; ++idx) {
                tempTruthProbs[idx] = sampleVariantTruth[idx] ? 1F : 0F;
            }
            bestPossibleLoss = getLoss(tempTruthProbs, variantIndices).toDouble();
            System.out.println("Best possible " + name + " loss = " + bestPossibleLoss);
        }

        private int[] getMiniBatchSplits(final int[] variantIndices, final int maxMiniBatchSize, final int numFilterableSampleVariants) {
            if(maxMiniBatchSize < getNumSamples()) {
                throw new IllegalArgumentException("in " + name + " DataSubset: maxMiniBatchSize (" + maxMiniBatchSize +
                                                   ") is less than numSamples (" + getNumSamples() + ")");
            }
            if(numFilterableSampleVariants <= maxMiniBatchSize) {
                return null;
            }
            final int[] numFilterableSamplesPerVariant = Arrays.stream(variantIndices)
                .mapToLong(XGBoostMinGqVariantFilter.this::getNumTrainableSamples)
                .mapToInt(i -> (int) i)
                .toArray();
            final List<Integer> miniBatchSplits = new ArrayList<>(2 * numFilterableSampleVariants / maxMiniBatchSize);
            int runningCount = 0;
            for(int i = 0; i < variantIndices.length; ++i) {
                runningCount += numFilterableSamplesPerVariant[i];
                if(runningCount > maxMiniBatchSize) {
                    miniBatchSplits.add(i);
                    runningCount = numFilterableSamplesPerVariant[i];
                }
            }
            miniBatchSplits.add(variantIndices.length);
            return miniBatchSplits.stream().mapToInt(Integer::intValue).toArray();
        }

        public int size() { return numFilterableSampleVariants; }
        public boolean usesMiniBatches() { return miniBatchSplits != null; }
        public int getNumMiniBatches() { return miniBatchSplits == null ? 1 : miniBatchSplits.length; }

        private int[] getBatchVariantIndices(final int miniBatchIndex) {
            final int i1 = miniBatchIndex == 0 ? 0 : miniBatchSplits[miniBatchIndex - 1];
            final int i2 = miniBatchSplits[miniBatchIndex];
            return Arrays.stream(variantIndices, i1, i2).toArray();
        }

        private DMatrix getBatchDMatrix(final int miniBatchIndex) {
            return dMatrix == null ?
                getDMatrix(getBatchVariantIndices(miniBatchIndex)):
                dMatrix;
        }

        private float[][] getRawPredictions(final Booster booster, final DMatrix dMatrix) {
            try {
                return booster.predict(dMatrix, true, 0);
            } catch(XGBoostError xgBoostError) {
                throw new GATKException("In " + name + " DataSubset: Predict error", xgBoostError);
            }
        }

        private void displayFloatHistogram(final float[] arr, final String description, final int numBins) {
            double minValue = arr[0];
            double maxValue = arr[0];
            DoubleStream.Builder builder = DoubleStream.builder();
            for(final float val : arr) {
                if(val < minValue) {
                    minValue = val;
                } else if(val > maxValue) {
                    maxValue = val;
                }
                builder.add(val);
            }
            if(maxValue - minValue < 1e-3) {
                System.out.println(description);
                System.out.format("\t%f: 100%%\n", minValue);
            } else {
                displayHistogram(description, builder.build(), numBins, minValue, maxValue);
            }
        }

        private Loss getMiniBatchLoss(final Booster booster, final int miniBatchIndex) {
            final DMatrix miniBatchDMatrix = getBatchDMatrix(miniBatchIndex);
            final float[][] rawPredictions = getRawPredictions(booster, miniBatchDMatrix);
            final float[] pVSGood = pSampleVariantGood == null ? new float[rawPredictions.length] : pSampleVariantGood;
            for(int idx = 0; idx < rawPredictions.length; ++idx) {
                pVSGood[idx] = sigma(rawPredictions[idx][0]);
            }
            final Loss loss = getLoss(pVSGood, variantIndices);

            if(isTrainingSet) {
                calculateDerivatives(pVSGood);
                if(miniBatchIndex == 0) {
                    System.out.format("\tBoosting round %d on %s\n", 1 + getRound(), name);
                }
                try {
                    booster.boost(miniBatchDMatrix, lossDerivatives[0], lossDerivatives[1]);
                } catch(XGBoostError xgBoostError) {
                    throw new GATKException("In " + name + " DataSubset: Boost error", xgBoostError);
                }
            }

            return loss;
        }

        public DoubleStream streamPredictions(final Booster booster) {
            final DoubleStream.Builder builder = DoubleStream.builder();
            final int numMiniBatches = getNumMiniBatches();
            for(int miniBatchIndex = 0; miniBatchIndex < numMiniBatches; ++miniBatchIndex) {
                final DMatrix miniBatchDMatrix = getBatchDMatrix(miniBatchIndex);
                final float[][] rawPredictions = getRawPredictions(booster, miniBatchDMatrix);
                for (float[] rawPrediction : rawPredictions) {
                    builder.add(sigma(rawPrediction[0]));
                }
            }
            return builder.build();
        }

        public void predictOneRound(final Booster booster) {
            final Loss roundLoss;
            if(miniBatchSplits == null) {
                roundLoss = getMiniBatchLoss(booster, 0).divide(this.size());
            } else {
                roundLoss = IntStream.range(0, miniBatchSplits.length)
                    .mapToObj(idx -> this.getMiniBatchLoss(booster, idx))
                    .reduce(Loss::add)
                    .orElseThrow(RuntimeException::new)
                    .divide(this.size());
            }
            appendScore(roundLoss.toFloat());
            if(printProgress > 0) {
                System.out.format("\t%s: %.3f\n",name, roundLoss.toFloat());
            }
        }

        void calculateDerivatives(final float[] pGood) {
            final float[] d1Loss, d2Loss;
            if (miniBatchSplits == null) {
                d1Loss = this.d1Loss;
                d2Loss = this.d2Loss;
            } else {
                d1Loss = new float[pGood.length];
                d2Loss = new float[pGood.length];
                this.lossDerivatives[0] = d1Loss;
                this.lossDerivatives[1] = d2Loss;
            }
            float scale = 1F; // pGood.length;
            for (int i = 0; i < pGood.length; ++i) {
                final float pTrue = pGood[i];
                final float pFalse = 1F - pTrue;
                d1Loss[i] = sampleVariantTruth[i] ? -pFalse / scale : pTrue / scale;
                d2Loss[i] = pTrue * pFalse / scale;
            }
        }

        public float getLastScore() { return scores.get(scores.size() - 1); }

        public void appendScore(final float score) {
            scores.add(score);
            if(score < bestScore) {
                bestScore = score;
                bestScoreInd = getRound();
            }
        }

        public int getRound() {
            return scores.size() - 1;
        }

        public boolean isBestScore() {
            return bestScoreInd == getRound();
        }

        public boolean stop() {
            return getRound() == maxTrainingRounds || stopEarly();
        }

        public boolean stopEarly() {
            return earlyStoppingRounds > 0 && getRound() - bestScoreInd > earlyStoppingRounds;
        }
    }

    @Override
    protected boolean needsNormalizedProperties() { return false; }

    @Override
    protected float predict(float[] variantProperties) {
        try {
                return sigma(
                    booster.predict(getDMatrix(variantProperties), true, 0)[0][0]
                );
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error predicting", xgBoostError);
        }
    }

    @Override
    protected void loadModel(final InputStream inputStream) {
        try {
            booster = XGBoost.loadModel(inputStream);
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error loading XGBoost model", xgBoostError);
        }
    }

    @Override
    protected void saveModel(final OutputStream outputStream) {
        try {
            booster.saveModel(outputStream);
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error saving XGBoost model", xgBoostError);
        }
    }

    private Booster initializeBooster(final List<DataSubset> dataSubsets) {
        final DataSubset trainingSubset = dataSubsets.stream()
            .filter(dataSubset -> dataSubset.isTrainingSet)
            .findFirst()
            .orElseThrow(() -> new GATKException("dataSubsets does not contain a training set"));
        final Map<String, DMatrix> watches = new HashMap<>();
        // add watches for any DataSubset that a) is not the main training set, and
        //                                     b) does not use mini batches (it can comfortably hang out in memory)
        for(final DataSubset dataSubset : dataSubsets) {
            if(!(dataSubset.equals(trainingSubset) || dataSubset.usesMiniBatches())) {
                watches.put(dataSubset.name, dataSubset.dMatrix);
            }
        }

        try {
            // Do 0 rounds of training on the first mini-batch (if using, otherwise whole dMatrix) of training DataSubset
            return XGBoost.train(trainingSubset.getBatchDMatrix(0), getXgboostParams(), 0, watches, null, null);
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error creating Booster", xgBoostError);
        }
    }

    private boolean trainOneRound(final Booster booster, final List<DataSubset> dataSubsets) {
        // evaluate booster on all data sets, calculate derivatives on training set
        if(printProgress > 0) {
            System.out.println("Training round " + (dataSubsets.get(0).getRound() + 1));
        }

        boolean needSaveCheckpoint = !dataSubsets.isEmpty();
        boolean needStop = !dataSubsets.isEmpty();
        for(final DataSubset dataSubset : dataSubsets) {
            dataSubset.predictOneRound(booster);
            if(!dataSubset.isTrainingSet) {
                // check if need to stop iterating or save model checkpoint
                needSaveCheckpoint = needSaveCheckpoint && dataSubset.isBestScore();
                needStop = needStop && dataSubset.stop();
            }
        }

        // check if booster needs to be saved, or if early stopping is necessary
        if(needSaveCheckpoint) {
            saveModelCheckpoint();
        }
        return !needStop;
    }

    void displayFeatureImportance(final Booster booster) {
        final List<String> propertyNames = propertiesTable.getPropertyNames();
        final Map<String, Integer> featureScore;
        try {
            featureScore = booster.getFeatureScore((String) null).entrySet().stream()
                .collect(
                    Collectors.toMap(
                        entry -> propertyNames.get(Integer.parseInt(entry.getKey().substring(1))),
                        Map.Entry::getValue
                    )
                );
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error getting feature score", xgBoostError);
        }

        System.out.println("Feature importance:");
        featureScore.entrySet().stream()
            .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
            .forEachOrdered(entry -> System.out.format("\t%s: %d\n", entry.getKey(), entry.getValue()));
    }

    @Override
    protected void trainFilter() {
        final List<DataSubset> dataSubsets = new ArrayList<DataSubset>() {
            private static final long serialVersionUID = 0L;
            {
                add(new DataSubset(TRAIN_MAT_KEY, getTrainingIndices(), true, maxMiniBatchSize));
                add(new DataSubset(VALIDATION_MAT_KEY, getValidationIndices(), false, maxMiniBatchSize));
            }
        };
        booster = initializeBooster(dataSubsets);
        //noinspection StatementWithEmptyBody
        while (trainOneRound(booster, dataSubsets))
            ;

        loadModelCheckpoint();
        displayHistogram("Final training adjusted GQ histogram",
                          dataSubsets.get(0).streamPredictions(booster).mapToInt(p -> phredScale(1.0 - p)),true);
        displayHistogram("Final training probability histogram",
                          dataSubsets.get(0).streamPredictions(booster),20, 0.0, 1.0);
        displayFeatureImportance(booster);
    }
}
