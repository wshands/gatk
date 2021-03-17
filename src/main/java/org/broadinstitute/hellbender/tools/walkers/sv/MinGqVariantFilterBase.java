package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.JSONValue;
import net.minidev.json.parser.ParseException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.samples.PedigreeValidationType;
import org.broadinstitute.hellbender.utils.samples.SampleDBBuilder;
import org.broadinstitute.hellbender.utils.samples.Trio;

import java.io.*;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.*;

import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.jetbrains.annotations.NotNull;

import static org.apache.commons.math3.util.FastMath.*;


/**
 * Extract matrix of properties for each variant.
 * Also extract, numVariants x numTrios x 3 tensors of allele count and genotype quality.
 * These data will be used to train a variant filter based on min GQ" (and stratified by other variant properties) that
 * maximizes the admission of variants with Mendelian inheritance pattern while omitting non-Mendelian variants.
 *
 * Derived class must implement abstract method trainFilter()
 */
public abstract class MinGqVariantFilterBase extends VariantWalker {
    @Argument(fullName=StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName=StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME,
            doc="Pedigree file, necessary for \"TRAIN\" mode, ignored in \"FILTER\" mode.", optional = true)
    public GATKPath pedigreeFile = null;

    @Argument(fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc="Output file. In \"FILTER\" mode it is the filtered VCF. In \"TRAIN\" mode it is the saved final filtering model.")
    public GATKPath outputFile;

    @Argument(fullName="model-file", shortName="m", optional=true,
              doc="Path to saved pre-existing filter model. In \"FILTER\" mode, this is a mandatory argument."
                 +" In \"TRAIN\" mode, this is an optional argument for additional training based on an existing model.")
    public GATKPath modelFile = null;

    @Argument(fullName="properties-table-file", optional=true,
            doc="Path to save properties table as JSON for analysis outside of this program. Only has an effect in TRAIN mode.")
    public GATKPath propertiesTableFile = null;

    @Argument(fullName="truth-file", shortName="t", optional=true,
              doc="Path to JSON file with truth data. Keys are sample IDs and values are objects with key \"good\""
                 +" corresponding to a list of known true variant IDs, and key \"bad\" corresponding to a list of known"
                 +" bad variant IDs")
    public GATKPath truthFile = null;

    private enum RunMode { TRAIN, FILTER }
    @Argument(fullName="mode", doc="Mode of operation: either \"TRAIN\" or \"FILTER\"")
    public RunMode runMode;

    @Argument(fullName="keep-homvar", shortName="kh", doc="Keep homvar variants even if their GQ is less than min-GQ", optional=true)
    public boolean keepHomvar = true;

    @Argument(fullName="keep-multiallelic", shortName="km", doc="Keep multiallelic variants even if their GQ is less than min-GQ", optional=true)
    public boolean keepMultiallelic = true;

    @Argument(fullName="validation-proportion", shortName="vp", doc="Proportion of variants to set aside for cross-validation",
              optional=true, minValue=0.0, maxValue=1.0)
    public double validationProportion = 0.5;

    @Argument(fullName="report-min-gq-filter-threshold", shortName="rf", optional=true, minValue=0.0, maxValue=1.0,
              doc="Add \"" + EXCESSIVE_MIN_GQ_FILTER_KEY + "\" to FILTER if the proportion of samples with calls filtered by "
                 + MIN_GQ_KEY + " is greater than this threshold.")
    public double reportMinGqFilterThreshold = 0.005;

    @Argument(fullName="truth-weight", shortName="tw", optional=true, minValue=0.0, maxValue=1.0,
            doc="Weight for truth data in loss function.")
    public static double truthWeight = 0.5;

    @Argument(fullName="print-progress", shortName="p", doc="Print progress of fit during training", optional = true)
    public int printProgress = 1;

    static final String minSamplesToEstimateAlleleFrequencyKey = "min-samples-to-estimate-allele-frequency";
    @Argument(fullName=minSamplesToEstimateAlleleFrequencyKey, shortName="ms", optional=true,
              doc="If the VCF does not have allele frequency, estimate it from the sample population if there are at least this many samples. Otherwise throw an exception.")
    public int minSamplesToEstimateAlleleFrequency = 100;

    @Argument(fullName="genome-tract", shortName="gt", optional=true)
    final List<String> genomeTractFiles = new ArrayList<>();

    @Argument(fullName="minGQ", shortName="gq", optional=true)
    int minGQ = FastMath.round(phredScale(0.5));

    @Argument(fullName="prediction-scale-factor", doc="Scale factor for raw predictions to logits", optional=true)
    public double predictionScaleFactor = 1.0;


    List<TractOverlapDetector> tractOverlapDetectors = null;

    static final Map<String, double[]> propertyBinsMap = new HashMap<String, double[]>() {
        private static final long serialVersionUID = 0;
        {
            put(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {0.05, 0.5});
            put(SVLEN_KEY, new double[] {2.0, 500.0, 5000.0});
        }
    };



    protected List<String> trainingSampleIds = null; // numSamples list of IDs for samples that will be used for training
    protected int[][] trioSampleIndices = null; // numTrios x 3 matrix of sample indices (paternal, maternal, child) for each trio
    protected Set<Integer> allTrioSampleIndices; // set of all sample indices that are in a trio
    protected final Map<String, Integer> trainingSampleIndices = new HashMap<>(); // map from sampleId to numerical index, for samples used in training
    protected final PropertiesTable propertiesTable = new PropertiesTable();


    // numVariants x numSamples matrix of allele counts:
    protected List<int[]> sampleVariantAlleleCounts = new ArrayList<>();
    // numVariants x numSamples matrix of genotype qualities:
    protected List<int[]> sampleVariantGenotypeQualities = new ArrayList<>();
    // map from variantIndex to array of known good GQ values for that variant
    protected Map<Integer, int[]> goodVariantGqs = null;
    // map from variantIndex to array of known good/bad sample indices for that variant
    protected final Map<Integer, Set<Integer>> goodSampleVariantIndices = new HashMap<>();
    protected final Map<Integer, Set<Integer>> badSampleVariantIndices = new HashMap<>();
    // map from variantIndex to array of known bad GQ values for that variant
    protected Map<Integer, int[]> badVariantGqs = null;

    // numVariants array of property-category bins
    protected int[] propertyBins = null;
    int numPropertyBins;
    protected long[] variantsPerPropertyBin = null;
    protected String[] propertyBinLabels = null;
    protected double[] propertyBinInheritanceWeights = null;
    protected double[] propertyBinTruthWeights = null;

    protected Random randomGenerator = Utils.getRandomGenerator();

    private VariantContextWriter vcfWriter = null;

    private int numVariants;
    private int numSamples;
    private int numTrios;

    private static final float PROB_EPS = 1.0e-3F;
    private static final Set<String> GAIN_SV_TYPES = new HashSet<>(Arrays.asList("INS", "DUP", "MEI", "ITX"));
    private static final Set<String> BREAKPOINT_SV_TYPES = new HashSet<>(Arrays.asList("BND", "CTX"));
    private static final double SV_EXPAND_RATIO = 1.0; // ratio of how much to expand SV gain range to SVLEN
    private static final int BREAKPOINT_HALF_WIDTH = 50; // how wide to make breakpoint intervals for tract overlap
    private static final String[] DISTAL_TARGETS_KEYS = {"SOURCE", "CPX_INTERVALS"}; // keys for distal targets in VCF entries
    private static final String SVLEN_KEY = "SVLEN";
    private static final String EVIDENCE_KEY = "EVIDENCE";
    private static final String FILTER_KEY = "FILTER";
    private static final String NO_EVIDENCE = "NO_EVIDENCE";
    private static final String ALGORITHMS_KEY = "ALGORITHMS";
    private static final String NO_ALGORITHM = "NO_ALGORITHM";
    private static final String MIN_GQ_KEY = "MINGQ";
    private static final String EXCESSIVE_MIN_GQ_FILTER_KEY = "LOW_QUALITY";
    private static final String MULTIALLELIC_FILTER = "MULTIALLELIC";
    private static final String GOOD_VARIANT_TRUTH_KEY = "good";
    private static final String BAD_VARIANT_TRUTH_KEY = "bad";
    private static final String CHR2_KEY = "CHR2";
    private static final String POS2_KEY = "POS2";
    private static final String END2_KEY = "END2";

    private static final String PE_GQ_KEY = "PE_GQ";
    private static final String RD_GQ_KEY = "RD_GQ";
    private static final String SR_GQ_KEY = "SR_GQ";
    private static final int MISSING_GQ_VAL = -1;

    // properties used to gather main matrix / tensors during apply()
    private Set<Trio> pedTrios = null;
    private Map<String, Set<String>> goodSampleVariants = null;
    private Map<String, Set<String>> badSampleVariants = null;

    // train/validation split indices
    private int[] trainingIndices;
    private int[] validationIndices;

    // stats on tool actions
    private int numFilteredGenotypes;

    // properties for calculating f1 score or estimating its pseudo-derivatives
    private Boolean [][] isTrueVariant; // numVariants x numSamples
    private int[] perVariantOptimalMinGq = null; // <- REMOVE
    protected int[] maxDiscoverableMendelianAc = null;
    long numDiscoverableMendelianAc = 0;
    protected double optimalProportionOfSampleVariantsPassing;

    protected final int getNumVariants() { return propertiesTable.getNumRows(); }
    protected final int getNumSamples() { return numSamples; }
    protected final int getNumTrios() { return numTrios; }
    protected final int getNumProperties() { return propertiesTable.getNumProperties(); }
    protected final int[] getTrainingIndices() { return trainingIndices; }
    protected final int[] getValidationIndices() { return validationIndices; }

    /**
     * Entry-point function to initialize the samples database from input data
     */
    private void getPedTrios() {
        if(pedigreeFile == null) {
            throw new UserException.BadInput(StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME + " must be specified in \"TRAIN\" mode");
        }
        final SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        sampleDBBuilder.addSamplesFromPedigreeFiles(Collections.singletonList(pedigreeFile));

        final Set<String> vcfSamples = new HashSet<>(getHeaderForVariants().getGenotypeSamples());
        // get the trios from the pedigree file, keeping only those that are fully present in the input VCF
        pedTrios = sampleDBBuilder.getFinalSampleDB()
            .getTrios()
            .stream()
            .filter(trio -> vcfSamples.contains(trio.getPaternalID()) && vcfSamples.contains(trio.getMaternalID())
                            && vcfSamples.contains(trio.getChildID()))
            .collect(Collectors.toSet());
        if(pedTrios.isEmpty()) {
            throw new UserException.BadInput(
                "The pedigree file (" + pedigreeFile + ") does not contain any trios that are present in the input VCF ("
                + drivingVariantFile + ")"
            );
        }
        numTrios = pedTrios.size();
        // collect ped trios into training samples
        pedTrios.stream()
                .flatMap(trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID()))
                .forEach(this::addTrainingSampleId);
        // keep track of the samples that are in trios (those will always have "truth" through inheritance)
        trioSampleIndices = pedTrios.stream()
                .map(
                    trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID())
                            .mapToInt(trainingSampleIndices::get)
                            .toArray()
                )
                .toArray(int[][]::new);

        allTrioSampleIndices = pedTrios.stream()
                .flatMap(trio -> Stream.of(trio.getPaternalID(), trio.getMaternalID(), trio.getChildID()))
                .map(trainingSampleIndices::get)
                .collect(Collectors.toSet());

    }

    private void loadVariantTruthData() {
        if(truthFile == null) {
            System.out.println("No truth file specified. Not using truth data.");
            return;
        } else {
            System.out.println("Loading truth data from " + truthFile);
        }
        final JSONObject jsonObject;
        try (final InputStream fileInputStream = truthFile.getInputStream()){
            jsonObject = (JSONObject) JSONValue.parseWithException(fileInputStream);
        } catch (IOException | ParseException ioException) {
            throw new GATKException("Unable to parse JSON from inputStream", ioException);
        }
        final Set<String> vcfSamples = new HashSet<>(getHeaderForVariants().getGenotypeSamples());
        goodSampleVariants = new HashMap<>();
        badSampleVariants = new HashMap<>();
        for(final Map.Entry<String, Object> sampleTruthEntry : jsonObject.entrySet()) {
            final String sampleId = sampleTruthEntry.getKey();
            if(!vcfSamples.contains(sampleId)) {
                continue; // don't bother storing truth data that doesn't overlap the input VCF
            }
            final JSONObject sampleTruth = (JSONObject)sampleTruthEntry.getValue();
            for(final Object variantIdObj : (JSONArray)sampleTruth.get(GOOD_VARIANT_TRUTH_KEY)) {
                final String variantId = (String)variantIdObj;
                if(goodSampleVariants.containsKey(variantId)) {
                    goodSampleVariants.get(variantId).add(sampleId);
                } else {
                    goodSampleVariants.put(variantId, new HashSet<>(Collections.singleton(sampleId)) );
                }
            }
            for(final Object variantIdObj : (JSONArray)sampleTruth.get(BAD_VARIANT_TRUTH_KEY)) {
                final String variantId = (String)variantIdObj;
                if(badSampleVariants.containsKey(variantId)) {
                    badSampleVariants.get(variantId).add(sampleId);
                } else {
                    badSampleVariants.put(variantId, new HashSet<>(Collections.singleton(sampleId)));
                }
            }
        }
        if(goodSampleVariants.isEmpty() && badSampleVariants.isEmpty()) {
            System.out.println("Truth file specified (" + truthFile + "), but no samples/variants overlap with input VCF ("
                               + drivingVariantFile + "). Not using truth data.");
            goodSampleVariants = null;
            badSampleVariants = null;
            return;
        }
        // Add any new samples that have truth but are not in trios file to the training samples
        goodSampleVariants.values().stream().flatMap(Collection::stream).forEach(this::addTrainingSampleId);
        badSampleVariants.values().stream().flatMap(Collection::stream).forEach(this::addTrainingSampleId);
        // prepare to hold data for scoring
        goodVariantGqs = new HashMap<>();
        badVariantGqs = new HashMap<>();
    }

    private void addTrainingSampleId(final String sampleId) {
        trainingSampleIndices.put(sampleId, trainingSampleIndices.getOrDefault(sampleId, trainingSampleIndices.size()));
    }

    private void setTrainingSampleIds() {
        trainingSampleIds = trainingSampleIndices.entrySet()
            .stream()
            .sorted(Comparator.comparingInt(Map.Entry::getValue))
            .map(Map.Entry::getKey)
            .collect(Collectors.toList());
    }

    private void initializeVcfWriter() {
        vcfWriter = createVCFWriter(outputFile);
        final Set<VCFHeaderLine> hInfo = new LinkedHashSet<>(getHeaderForVariants().getMetaDataInInputOrder());
        final String filterableVariant = (keepMultiallelic ? "biallelic " : "") + (keepHomvar ? "het " : "") + "variant";
        hInfo.add(new VCFInfoHeaderLine(MIN_GQ_KEY, 1, VCFHeaderLineType.Integer,
                             "Minimum passing GQ for each " + filterableVariant));
        hInfo.add(new VCFFilterHeaderLine(EXCESSIVE_MIN_GQ_FILTER_KEY,
                               "More than " + (100 * reportMinGqFilterThreshold) + "% of sample GTs were masked as no-call GTs due to low GQ"));
        vcfWriter.writeHeader(new VCFHeader(hInfo, getHeaderForVariants().getGenotypeSamples()));
    }

    @Override
    public void onTraversalStart() {
        loadTrainedModel();  // load model and saved properties stats
            tractOverlapDetectors = genomeTractFiles.stream()
                .map(genomeTractFile -> {
                        try {
                            return new TractOverlapDetector(genomeTractFile);
                        } catch(IOException ioException) {
                            throw new GATKException("Error loading " + genomeTractFile, ioException);
                        }
                    })
                .collect(Collectors.toList());
        if(runMode == RunMode.TRAIN) {
            getPedTrios();  // get trios from pedigree file
            loadVariantTruthData(); // load variant truth data from JSON file
            setTrainingSampleIds(); // set data organizing training samples
            numSamples = trainingSampleIndices.size();
        } else {
            initializeVcfWriter();  // initialize vcfWriter and write header
            numFilteredGenotypes = 0;
            propertiesTable.setNumAllocatedRows(1);
        }
    }

    private static boolean mapContainsTrio(final Map<String, Integer> map, final Trio trio) {
        return map.containsKey(trio.getPaternalID()) && map.containsKey(trio.getMaternalID())
                && map.containsKey(trio.getChildID());
    }

    private static int[] getMappedTrioProperties(final Map<String, Integer> map, final Trio trio) {
        return new int[] {map.get(trio.getPaternalID()), map.get(trio.getMaternalID()), map.get(trio.getChildID())};
    }

    private boolean getIsMultiallelic(final VariantContext variantContext) {
        return variantContext.getNAlleles() > 2 || variantContext.getFilters().contains(MULTIALLELIC_FILTER);
    }

    private boolean getVariantIsFilterable(final VariantContext variantContext, final Map<String, Integer> sampleAlleleCounts) {
        final boolean maybeFilterable = !(keepMultiallelic && getIsMultiallelic(variantContext));
        if(maybeFilterable) {
            if(runMode == RunMode.FILTER) {
                // filter no matter what, because end-user may be interested in minGQ
                return true;
            } else {
                // check if any of the allele counts can be filtered so as not to train on unfilterable variants
                return pedTrios.stream()
                        .filter(trio -> mapContainsTrio(sampleAlleleCounts, trio))
                        .flatMapToInt(trio -> Arrays.stream(getMappedTrioProperties(sampleAlleleCounts, trio)))
                        .anyMatch(keepHomvar ? ac -> ac == 1 : ac -> ac > 0);
            }
        } else {
            return false;
        }
    }

    /**
     * form numTrios x 3 matrix of allele counts for specified variantIndex
     * @param variantIndex
     * @return
     */
    protected int[][] getTrioAlleleCountsMatrix(final int variantIndex) {
        final int[] sampleAlleleCounts = sampleVariantAlleleCounts.get(variantIndex);
        return Arrays.stream(trioSampleIndices)
            .map(
                trioIndices -> Arrays.stream(trioIndices)
                                .map(idx -> sampleAlleleCounts[idx])
                                .toArray()
            )
            .toArray(int[][]::new);
    }

    /**
     * form numTrios x 3 matrix of genotype qualities for specified variantIndex
     * @param variantIndex
     * @return
     */
    protected int[][] getTrioGenotypeQualitiesMatrix(final int variantIndex) {
        final int[] sampleGenotypeQualities = sampleVariantGenotypeQualities.get(variantIndex);
        return Arrays.stream(trioSampleIndices)
                .map(
                        trioIndices -> Arrays.stream(trioIndices)
                                .map(idx -> sampleGenotypeQualities[idx])
                                .toArray()
                )
                .toArray(int[][]::new);
    }


    protected boolean getSampleVariantIsFilterable(final int variantIndex, final int sampleIndex) {
        return alleleCountIsFilterable(sampleVariantAlleleCounts.get(variantIndex)[sampleIndex]);
    }

    /**
     * A sample variant is filterable if it's het, or if it's HOMVAR and keepHomVar is not true
     */
    protected boolean alleleCountIsFilterable(final int alleleCount) {
        return keepHomvar ? alleleCount == 1 : alleleCount > 0;
    }

    /**
     * A sample Variant is trainable if a) it is filterable, and b) there is truth data (inheritance or known good/bad)
     */
    protected boolean getSampleVariantIsTrainable(final int variantIndex, final int sampleIndex) {
        if(allTrioSampleIndices.contains(sampleIndex) ||
           goodSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet()).contains(sampleIndex) ||
           badSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet()).contains(sampleIndex)) {
            return getSampleVariantIsFilterable(variantIndex, sampleIndex);
        } else {
            return false;
        }
    }

    private boolean hasDistalTarget(final VariantContext variantContext) {
        return variantContext.hasAttribute(CHR2_KEY) && variantContext.hasAttribute(END2_KEY) && variantContext.hasAttribute(POS2_KEY);
    }

    private SimpleInterval getDistalTarget(final String targetString) {
        return new SimpleInterval(targetString.substring(targetString.indexOf("_") + 1));
    }

    private void getTractProperties(final TractOverlapDetector tractOverlapDetector,
                                    final VariantContext variantContext) {
        // get range of variant (use expanded location for insertions or duplications)

//        To-do:
//            1) Change this to form list of locations
//                -> list of overlaps -> min / max
//                -> spans by looping over all possible pairs if appropriate
//                list of locations
//                    BND-type: left + right
//                    all others have center
//                    if end - pos > 1 then add left + right
//                    distal targets get added in
//            2) compare to score with "correct" answers from initial iteration
//                weight by bin
//                weight by gq - optimal min gq
        final List<Locatable> overlapCheckLocations = new ArrayList<>();
        final String svType = variantContext.getAttributeAsString(VCFConstants.SVTYPE, null);
        final int centerInd;
        if(BREAKPOINT_SV_TYPES.contains(svType)) {
            // Add checks at the breakpoint locations. They are points so we want to check tract overlaps in the regions
            // around them
            overlapCheckLocations.add(
                new SimpleInterval(variantContext.getContig(), variantContext.getStart() - BREAKPOINT_HALF_WIDTH,
                             variantContext.getStart() + BREAKPOINT_HALF_WIDTH)
            );
            overlapCheckLocations.add(
                new SimpleInterval(variantContext.getAttributeAsString(CHR2_KEY, null),
                             variantContext.getEnd() - BREAKPOINT_HALF_WIDTH,
                              variantContext.getEnd() + BREAKPOINT_HALF_WIDTH)
            );
        } else {
            // Insertions can have zero width but they have uncertainty in their true location. Expand their central
            // interval for the purpose of checking tract overlap. Do this for any type that's "narrow"
            final int width = variantContext.getEnd() - variantContext.getStart();
            final boolean is_narrow = width < 2 * BREAKPOINT_HALF_WIDTH;
            final int expand = is_narrow ?
                max(
                    (int)(abs(variantContext.getAttributeAsInt(SVLEN_KEY, 0)) * SV_EXPAND_RATIO),
                    BREAKPOINT_HALF_WIDTH
                ) :
                0;
            overlapCheckLocations.add(
                new SimpleInterval(variantContext.getContig(),
                    variantContext.getStart() - expand,
                    variantContext.getStart() + expand)
            );
            // if the interval isn't narrow, add the end points as potential breakpoints to check
            if(!is_narrow) {
                overlapCheckLocations.add(
                    new SimpleInterval(variantContext.getContig(),
                                      variantContext.getStart() - BREAKPOINT_HALF_WIDTH,
                                      variantContext.getStart() + BREAKPOINT_HALF_WIDTH)
                );
                overlapCheckLocations.add(
                    new SimpleInterval(variantContext.getContig(),
                                 variantContext.getEnd() - BREAKPOINT_HALF_WIDTH,
                                 variantContext.getEnd() + BREAKPOINT_HALF_WIDTH)
                );
            }
        }
        // Check for other distal locations in "SOURCE" field
        //   they may be null (filter out)
        //   they may be of form SVTYPE_SIMPLEINTERVAL (chop off everything before "_")
        Arrays.stream(DISTAL_TARGETS_KEYS)
                .map(key -> variantContext.getAttributeAsStringList(key, null))
                .filter(Objects::nonNull)
                .flatMap(List::stream)
                .map(targetString -> new SimpleInterval(targetString.substring(targetString.indexOf("_") + 1)))
                .forEach(overlapCheckLocations::add);

        final double[] overlaps;
        if(tractOverlapDetector.hasOther()) { // This tract has paired intervals
            overlaps = overlapCheckLocations.stream()
                    .mapToDouble(location -> tractOverlapDetector.getPrimaryOverlapFraction(location)
                                             + tractOverlapDetector.getOtherOverlapfraction(location))
                    .toArray();
            // check for spanning by streaming over all possible pairs
            final boolean spans = IntStream.range(0, overlapCheckLocations.size() - 1)
                .anyMatch(
                    i -> IntStream.range(i + 1, overlapCheckLocations.size()).anyMatch(
                        j -> tractOverlapDetector.spansPrimaryAndOther(overlapCheckLocations.get(i),
                                                                       overlapCheckLocations.get(j))
                    )
                );
            propertiesTable.append(
                tractOverlapDetector.getName() + "_spans",
                spans
            );
        } else {  // this tract has simple intervals
            overlaps = overlapCheckLocations.stream()
                .mapToDouble(tractOverlapDetector::getPrimaryOverlapFraction)
                .toArray();
        }
        // breakpoint types have no "center", otherwise get overlap of primary location
        propertiesTable.append(
            tractOverlapDetector.getName() + "_center",
            BREAKPOINT_SV_TYPES.contains(svType) ? 0.0 : overlaps[0]
        );
        propertiesTable.append(
            tractOverlapDetector.getName() + "_min",
            Arrays.stream(overlaps).min().orElse(0.0)
        );
        propertiesTable.append(
            tractOverlapDetector.getName() + "_max",
            Arrays.stream(overlaps).max().orElse(0.0)
        );
    }

    public int getGenotypeAttributeAsInt(final Genotype genotype, final String key, Integer defaultValue) {
        Object x = genotype.getExtendedAttribute(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) {
            if(defaultValue == null) {
                throw new IllegalArgumentException("Genotype is missing value of " + key);
            } else {
                return defaultValue;
            }
        }
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.parseInt((String)x); // throws an exception if this isn't a string
    }

    private int[] getGenotypeAttributeAsInt(final Iterable<Genotype> sampleGenotypes, final String attributeKey,
                                            final int missingAttributeValue) {
        return StreamSupport.stream(sampleGenotypes.spliterator(), false)
                        .mapToInt(g -> getGenotypeAttributeAsInt(g, attributeKey, missingAttributeValue))
                        .toArray();
    }

    /**
     * Accumulate properties for variant matrix, and allele counts, genotype quality for trio tensors
     */
    @Override
    public void apply(VariantContext variantContext, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        // get per-sample allele counts as a map indexed by sample ID
        final Map<String, Integer> sampleAlleleCounts = variantContext.getGenotypes().stream().collect(
                Collectors.toMap(
                        Genotype::getSampleName,
                        g -> g.getAlleles().stream().mapToInt(a -> a.isReference() ? 0 : 1).sum()
                )
        );
        if(!getVariantIsFilterable(variantContext, sampleAlleleCounts)) {
            // no need to train on unfilterable variants, and filtering is trivial
            if(runMode == RunMode.FILTER) {
                ++numVariants;
                vcfWriter.add(variantContext); // add variantContext unchanged
            }
            return;
        }
        ++numVariants;

        double alleleFrequency = variantContext.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
        if(alleleFrequency <= 0) {
            if(variantContext.getNSamples() <= minSamplesToEstimateAlleleFrequency) {
                throw new GATKException("VCF does not have " + VCFConstants.ALLELE_FREQUENCY_KEY + " annotated or enough samples to estimate it ("
                                        + minSamplesToEstimateAlleleFrequencyKey + "=" + minSamplesToEstimateAlleleFrequency + " but there are "
                                        + variantContext.getNSamples() + " samples)");
            }
            // VCF not annotated with allele frequency, guess it from allele counts
            final int numAlleles = variantContext.getGenotypes().stream().mapToInt(Genotype::getPloidy).sum();
            alleleFrequency = sampleAlleleCounts.values().stream().mapToInt(Integer::intValue).sum() / (double) numAlleles;
        }
        propertiesTable.append(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFrequency);

        final String svType = variantContext.getAttributeAsString(VCFConstants.SVTYPE, null);
        if(svType == null) {
            throw new GATKException("Missing " + VCFConstants.SVTYPE + " for variant " + variantContext.getID());
        }
        propertiesTable.append(VCFConstants.SVTYPE, svType);

        final int svLen = variantContext.getAttributeAsInt(SVLEN_KEY, Integer.MIN_VALUE);
        if(svLen == Integer.MIN_VALUE) {
            throw new GATKException("Missing " + SVLEN_KEY + " for variant " + variantContext.getID());
        }
        propertiesTable.append(SVLEN_KEY, svLen);

        propertiesTable.append(FILTER_KEY, variantContext.getFilters());

        final Set<String> vcEvidence = Arrays.stream(
                    variantContext.getAttributeAsString(EVIDENCE_KEY, NO_EVIDENCE)
                        .replaceAll("[\\[\\] ]", "").split(",")
            ).map(ev -> ev.equals(".") ? NO_EVIDENCE : ev).collect(Collectors.toSet());
        if(vcEvidence.isEmpty()) {
            throw new GATKException("Missing " + EVIDENCE_KEY + " for variant " + variantContext.getID());
        }
        propertiesTable.append(EVIDENCE_KEY, vcEvidence);


        final Set<String> vcAlgorithms = Arrays.stream(
                variantContext.getAttributeAsString(ALGORITHMS_KEY, NO_ALGORITHM)
                        .replaceAll("[\\[\\] ]", "").split(",")
        ).map(ev -> ev.equals(".") ? NO_ALGORITHM : ev).collect(Collectors.toSet());
        if(vcAlgorithms.isEmpty()) {
            throw new GATKException("Missing " + ALGORITHMS_KEY + " for variant " + variantContext.getID());
        }
        propertiesTable.append(ALGORITHMS_KEY, vcAlgorithms);

        for(final TractOverlapDetector tractOverlapDetector : tractOverlapDetectors) {
            getTractProperties(tractOverlapDetector, variantContext);
        }

        final Iterable<Genotype> sampleGenotypes = variantContext.getGenotypesOrderedBy(trainingSampleIds);
        propertiesTable.append(
            VCFConstants.GENOTYPE_QUALITY_KEY,
            StreamSupport.stream(sampleGenotypes.spliterator(), false).mapToInt(Genotype::getGQ).toArray()
        );

        Stream.of(PE_GQ_KEY, RD_GQ_KEY, SR_GQ_KEY).forEach(
            intGenotypeAttributeKey -> propertiesTable.append(
                    intGenotypeAttributeKey,
                    getGenotypeAttributeAsInt(sampleGenotypes, intGenotypeAttributeKey, MISSING_GQ_VAL)
            )
        );

        if(runMode == RunMode.TRAIN) {
            // get per-sample genotype qualities as a map indexed by sample ID
            final Map<String, Integer> sampleGenotypeQualities = variantContext.getGenotypes().stream().collect(
                    Collectors.toMap(Genotype::getSampleName, Genotype::getGQ)
            );

            if(goodSampleVariants != null) {
                if(goodSampleVariants.containsKey(variantContext.getID())) {
                    // Get sample IDs of known good variants that are filterable
                    final String[] knownGoodSampleIds = goodSampleVariants.get(variantContext.getID())
                            .stream()
                            .filter(keepHomvar ?
                                    sampleId -> sampleAlleleCounts.get(sampleId) == 1 :
                                    sampleId -> sampleAlleleCounts.get(sampleId) > 0)
                            .toArray(String[]::new);
                    if(knownGoodSampleIds.length > 0) {
                        // get indices and GQ values of good samples for this variant
                        goodSampleVariantIndices.put(
                            numVariants - 1, Arrays.stream(knownGoodSampleIds).map(trainingSampleIndices::get).collect(Collectors.toSet())
                        );
                        // Get GQ values of known good variants that are filterable
                        goodVariantGqs.put(
                            numVariants - 1, Arrays.stream(knownGoodSampleIds).mapToInt(sampleGenotypeQualities::get).toArray()
                        );
                    }
                }
                if(badSampleVariants.containsKey(variantContext.getID())) {
                    // Get sample IDs of known bad variants that are filterable
                    final String[] knownBadSampleIds = badSampleVariants.get(variantContext.getID())
                            .stream()
                            .filter(keepHomvar ?
                                    sampleId -> sampleAlleleCounts.get(sampleId) == 1 :
                                    sampleId -> sampleAlleleCounts.get(sampleId) > 0)
                            .toArray(String[]::new);
                    if(knownBadSampleIds.length > 0) {
                        // get indices and GQ values of good samples for this variant
                        badSampleVariantIndices.put(
                                numVariants - 1, Arrays.stream(knownBadSampleIds).map(trainingSampleIndices::get).collect(Collectors.toSet())
                        );
                        // Get GQ values of known good variants that are filterable
                        badVariantGqs.put(
                                numVariants - 1, Arrays.stream(knownBadSampleIds).mapToInt(sampleGenotypeQualities::get).toArray()
                        );
                    }
                }
            }

            // get allele counts for training samples for this variant
            sampleVariantAlleleCounts.add(
                trainingSampleIds.stream().mapToInt(sampleAlleleCounts::get).toArray()
            );

            // get genotype qualities for training samples for this variant
            sampleVariantGenotypeQualities.add(
                trainingSampleIds.stream().mapToInt(sampleGenotypeQualities::get).toArray()
            );
        } else {
            propertiesTable.validateAndFinalize();
            vcfWriter.add(filterVariantContext(variantContext));
            propertiesTable.clearRows();
        }
    }

    private VariantContext filterVariantContext(final VariantContext variantContext) {
        final Genotype[] genotypes = new Genotype[variantContext.getNSamples()];
        int numFiltered = 0;
        int sampleIndex = 0;
        for(final Genotype genotype : variantContext.getGenotypes()) {
            if(getSampleVariantIsFilterable(0, sampleIndex)) {
                final int gq = adjustedGq(
                    propertiesTable.getPropertiesRow(0, sampleIndex, needsNormalizedProperties())
                );
                if (gq >= minGQ) {
                    genotypes[sampleIndex] = new GenotypeBuilder(genotype).GQ(gq).make();

                } else {
                    genotypes[sampleIndex] = new GenotypeBuilder(genotype)
                        .alleles(GATKVariantContextUtils.noCallAlleles(genotype.getPloidy()))
                        .make();
                    ++numFiltered;
                }
            } else{
                genotypes[sampleIndex] = genotype;
            }
            ++sampleIndex;
        }
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variantContext)
            .genotypes(genotypes).attribute(MIN_GQ_KEY, minGQ);
        if(numFiltered > variantContext.getNSamples() * reportMinGqFilterThreshold) {
            variantContextBuilder.filter(EXCESSIVE_MIN_GQ_FILTER_KEY);
        }
        numFilteredGenotypes += numFiltered;

        return variantContextBuilder.make();
    }


    private boolean notCloseEnough(final String strRep, final double val, final double tol) {
        return FastMath.abs(Double.parseDouble(strRep) - val) > FastMath.abs(tol * val);
    }

    private String[] getPropertyBinNames( final PropertiesTable.Property property, final double[] bins) {
        double minVal = bins[0];
        double maxVal = bins[bins.length - 1];
        for(int variantIndex = 0; variantIndex < property.getNumRows(); ++variantIndex) {
            final float val = property.getAsFloat(variantIndex, 0, false);
            if(val < minVal) {
                minVal = val;
            } else if(val > maxVal) {
                maxVal = val;
            }
        }

        final String[] propertyBinNames = new String[bins.length + 1];
        for(int i = 0; i < bins.length + 1; ++i) {
            final double low = i == 0 ? minVal : bins[i - 1];
            final double high = i == bins.length ? maxVal : bins[i];
            int precision = 0;
            String lowStr = String.format("%." + precision + "f", low);
            String highStr = String.format("%." + precision + "f", high);
            while(lowStr.equals(highStr)) {
                ++precision;
                lowStr = String.format("%." + precision + "f", low);
                highStr = String.format("%." + precision + "f", high);
            }
            int lowPrecision = precision;
            while(notCloseEnough(lowStr, low, 0.1)) {
                ++lowPrecision;
                lowStr = String.format("%." + lowPrecision + "f", low);
            }
            int highPrecision = precision;
            while(notCloseEnough(highStr, high, 0.1)) {
                ++highPrecision;
                highStr = String.format("%." + highPrecision + "f", high);
            }
            propertyBinNames[i] = property.name + "=" + lowStr + "-" + highStr;
        }
        return propertyBinNames;
    }

    private void setPropertyBins() {
        propertyBins = new int[numVariants]; // guaranteed to be all zeros by Java spec

        final PropertiesTable.StringArrProperty svType = (PropertiesTable.StringArrProperty) propertiesTable.get(VCFConstants.SVTYPE);
        final List<String> allSvTypes = svType.getAllLabels();
        final Map<String, Integer> binNames = IntStream.range(0, allSvTypes.size())
            .boxed()
            .collect(Collectors.toMap(i -> VCFConstants.SVTYPE + "=" + allSvTypes.get(i), i -> i));
        IntStream.range(0, numVariants).forEach(variantIndex -> {
            propertyBins[variantIndex] = binNames.get(VCFConstants.SVTYPE + "=" + svType.getAsString(variantIndex));
        });

        propertyBinsMap.forEach((propertyName, bins) -> {
            final List<String> oldBinNames = binNames.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue())
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
            binNames.clear();

            final PropertiesTable.Property property = propertiesTable.get(propertyName);
            final String[] propertyBinNames = getPropertyBinNames(property, bins);
            IntStream.range(0, numVariants).forEach(variantIndex -> {
                int propertyBin = Arrays.binarySearch(bins, property.getAsFloat(variantIndex, 0, false));
                if(propertyBin < 0) { propertyBin = ~propertyBin; }
                final String newBinName = oldBinNames.isEmpty() ?
                    propertyBinNames[propertyBin] :
                    oldBinNames.get(propertyBins[variantIndex]) + ";" + propertyBinNames[propertyBin];
                if(binNames.containsKey(newBinName)) {
                    propertyBins[variantIndex] = binNames.get(newBinName);
                } else {
                    final int newBin = binNames.size();
                    binNames.put(newBinName, newBin);
                    propertyBins[variantIndex] = newBin;
                }
            });
        });
        numPropertyBins = binNames.size();
        propertyBinLabels = new String[numPropertyBins];
        binNames.forEach((key, value) -> propertyBinLabels[value] = key);
        variantsPerPropertyBin = new long[numPropertyBins];
        Arrays.stream(propertyBins).forEach(bin -> ++variantsPerPropertyBin[bin]);
    }

    private IntStream streamFilterableGq(final IntStream acStream, final IntStream gqStream) {
        final PrimitiveIterator.OfInt acIterator = acStream.iterator();
        return gqStream.filter(keepHomvar ? gc -> acIterator.nextInt() == 1 : gc -> acIterator.nextInt() > 0);
    }

    protected IntStream getCandidateMinGqs(final int variantIndex) {
        final IntStream alleleCountsStream = Arrays.stream(getTrioAlleleCountsMatrix(variantIndex)).flatMapToInt(Arrays::stream);
        final IntStream genotypeQualitiesStream = Arrays.stream(getTrioGenotypeQualitiesMatrix(variantIndex)).flatMapToInt(Arrays::stream);
        return getCandidateMinGqs(alleleCountsStream, genotypeQualitiesStream, null);
    }

    protected IntStream getCandidateMinGqs(final int[] rowIndices) {
        final IntStream alleleCountsStream = Arrays.stream(rowIndices).flatMap(
                rowIndex -> Arrays.stream(getTrioAlleleCountsMatrix(rowIndex)).flatMapToInt(Arrays::stream)
        );
        final IntStream genotypeQualitiesStream = Arrays.stream(rowIndices).flatMap(
                rowIndex -> Arrays.stream(getTrioGenotypeQualitiesMatrix(rowIndex)).flatMapToInt(Arrays::stream)
        );
        return getCandidateMinGqs(alleleCountsStream, genotypeQualitiesStream, null);
    }

    protected IntStream getCandidateMinGqs(final int[][] alleleCounts, final int[][] genotypeQualities,
                                           final Integer preexistingCandidateMinGq) {
        final IntStream alleleCountsStream = Arrays.stream(alleleCounts).flatMapToInt(Arrays::stream);
        final IntStream genotypeQualitiesStream = Arrays.stream(genotypeQualities).flatMapToInt(Arrays::stream);
        return getCandidateMinGqs(alleleCountsStream, genotypeQualitiesStream, preexistingCandidateMinGq);
    }

    protected IntStream getCandidateMinGqs(final IntStream alleleCountsStream, final IntStream genotypeQualitiesStream,
                                           final Integer preexistingCandidateMinGq) {
        // form list of candidate min GQ values that could alter filtering for this variant
        final List<Integer> candidateGq = streamFilterableGq(alleleCountsStream, genotypeQualitiesStream)
                .distinct()
                .sorted()
                .boxed()
                .collect(Collectors.toList());
        if(candidateGq.isEmpty()) {
            return null;
        }
        // consider filtering out all filterable variants by adding 1 to highest filterable GQ value
        candidateGq.add(1 + candidateGq.get(candidateGq.size() - 1));

        // find min GQ value that is redundant with preexistingCandidateMinGq
        final OptionalInt redundantIdx = preexistingCandidateMinGq == null ?
            OptionalInt.empty() :
            IntStream.range(0, candidateGq.size())
                    .filter(i -> candidateGq.get(i) >= preexistingCandidateMinGq).findFirst();

        // These candidate values are as large as they can be without altering filtering results.
        // If possible, move candidates down (midway to next change value) so that inaccuracy in predicting
        // candidateMinGq has tolerance in both directions.
        for(int index = candidateGq.size() - 2; index >= 0; --index) {
            final int gq = candidateGq.get(index);
            final int delta_gq = index == 0 ?
                    2:
                    gq - candidateGq.get(index - 1);
            if(delta_gq > 1) {
                candidateGq.set(index, gq - delta_gq / 2);
            }
        }

        // remove candidate gq that is redundant with preexistingCandidateMinGq
        if (redundantIdx.isPresent()) {
            candidateGq.remove(redundantIdx.getAsInt());
            if(candidateGq.isEmpty()) {
                return null;
            }
        }

        // return stream of primitive int
        return candidateGq.stream().mapToInt(Integer::intValue);
    }


    private boolean isMendelian(final int fatherAc, final int motherAc, final int childAc) {
        // child allele counts should not exhibit de-novo mutations nor be missing inherited homvar
        final int minAc = fatherAc / 2 + motherAc / 2;
        final int maxAc = (fatherAc > 0 ? 1 : 0) + (motherAc > 0 ? 1 : 0);
        return (minAc <= childAc) && (childAc <= maxAc);
    }

    private double getPMendelian(final int fatherAc, final int motherAc, final int childAc,
                                 final float pFather, final float pMother, final float pChild,
                                 final double[] d1P) {
        // Note:
        //  -inherit 0 or 1 alleles from each parent so pInherit0 = 1 - pInherit1
        //  -similarly dPInherit0 = -dPInherit1
        //  -diagonal of Hessian is always zero, so don't need d2P
        final double dPInherit1Father = fatherAc / 2.0;
        final double pInherit1Father = dPInherit1Father * pFather;
        final double pInherit0Father = 1.0 - pInherit1Father;
        final double dPInherit1Mother = motherAc / 2.0;
        final double pInherit1Mother = dPInherit1Mother * pMother;
        final double pInherit0Mother = 1.0 - pInherit1Mother;

        final double pRef = pInherit0Father * pInherit0Mother;
        switch(childAc) {
            case 0:
                d1P[0] = -dPInherit1Father * pInherit0Mother;
                d1P[1] = -dPInherit1Mother * pInherit0Father;
                d1P[2] = 0.0;  // pChild is irrelevant
                return pRef;
            case 1:
                final double pHet = pInherit1Father * pInherit0Mother + pInherit0Father * pInherit1Mother;
                d1P[0] = dPInherit1Father * (pChild * (1.0 - 2.0 * pInherit1Mother) - (1.0 - pChild) * pInherit0Mother);
                d1P[1] = dPInherit1Mother * (pChild * (1.0 - 2.0 * pInherit1Father) - (1.0 - pChild) * pInherit0Father);
                d1P[2] = pHet - pRef;
                return pChild * pHet + (1.0 - pChild) * pRef;
            case 2:
                final double pHom = pInherit1Father * pInherit1Mother;
                d1P[0] = dPInherit1Father * (pChild * pInherit1Mother - (1.0 - pChild) * pInherit0Mother);
                d1P[1] = dPInherit1Mother * (pChild * pInherit1Father - (1.0 - pChild) * pInherit0Father);
                d1P[2] = pHom - pRef;
                return pChild * pHom + (1.0 - pChild) * pRef;
            default:
                throw new IllegalArgumentException("Illegal allele count: " + childAc);
        }
    }

    final void setSampleIndexToPredictionIndex(final Integer[] sampleIndexToPredictionIndex,
                                               final int[] sampleAlleleCounts,
                                               final Set<Integer> goodSampleIndices,
                                               final Set<Integer> badSampleIndices,
                                               final int startPredictionIndex) {
        int predictionIndex = startPredictionIndex;
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            if(alleleCountIsFilterable(sampleAlleleCounts[sampleIndex]) &&
                    (allTrioSampleIndices.contains(sampleIndex) ||
                            goodSampleIndices.contains(sampleIndex) ||
                            badSampleIndices.contains(sampleIndex))
            ) {
                sampleIndexToPredictionIndex[sampleIndex] = predictionIndex;
                ++predictionIndex;
            } else {
                sampleIndexToPredictionIndex[sampleIndex] = null;
            }
        }
    }

    protected float logitsToProb(final float predict) {
        return 1.0F / (1.0F + (float)FastMath.exp(-predictionScaleFactor * predict));
    }

    protected float probToLogits(final float sigma) {
        return (float)(FastMath.log(sigma / (1F - sigma)) / predictionScaleFactor);
    }


    /**
     * Convert from dLoss / dProb to dLoss / dLogits
     * @param d1LossDProb
     * @param d2LossDProb
     * @param d1LossDPredict
     * @param d2LossDPredict
     * @param predictIndex
     */
    protected void setLossDerivs(final double prob, final double d1LossDProb, final double d2LossDProb,
                                 final double weight, final double[] d1LossDLogits, final double[] d2LossDLogits,
                                 final int sampleIndex) {
        final double d1ProbDLogit = prob * (1.0 - prob) * predictionScaleFactor;
        final double d2ProbDLogit = d1ProbDLogit * (1.0 - 2.0 * prob) * predictionScaleFactor;
        d1LossDLogits[sampleIndex] += weight * d1LossDProb * d1ProbDLogit;
        d2LossDLogits[sampleIndex] += weight * (d2LossDProb * d1ProbDLogit * d1ProbDLogit
                                                + d1LossDProb * d2ProbDLogit);
    }

    protected void setLossDerivs(final double prob, final double d1LossDProb, final double d2LossDProb,
                                 final double weight, final float[] d1LossDLogits, final float[] d2LossDLogits,
                                 final int sampleIndex) {
        final double d1ProbDLogit = prob * (1.0 - prob) * predictionScaleFactor;
        final double d2ProbDLogit = d1ProbDLogit * (1.0 - 2.0 * prob) * predictionScaleFactor;
        d1LossDLogits[sampleIndex] = (float)(d1LossDLogits[sampleIndex] + weight * d1LossDProb * d1ProbDLogit);
        d2LossDLogits[sampleIndex] = (float)(d2LossDLogits[sampleIndex]
                                             + weight * (d2LossDProb * d1ProbDLogit * d1ProbDLogit
                                                        + d1LossDProb * d2ProbDLogit));
    }

    private Loss getLoss(final float[] pSample, final double[] d1PSample, final double[] d2PSample,
                           final int startPredictionIndex, final Integer[] sampleIndexToPredictionIndex,
                           final double variantInheritanceWeight, final double variantTruthWeight,
                           final int variantIndex, final List<Integer> priorSampleIndices) {
        final int[] sampleAlleleCounts = sampleVariantAlleleCounts.get(variantIndex);
        final Set<Integer> goodSampleIndices = goodSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet());
        final Set<Integer> badSampleIndices = badSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet());
        setSampleIndexToPredictionIndex(sampleIndexToPredictionIndex, sampleAlleleCounts,
                                        goodSampleIndices, badSampleIndices, startPredictionIndex);

        final double inheritanceDerivWeight = variantInheritanceWeight * pSample.length / numTrios;
        final double[] d1PMendelianDPSample = new double[3];
        double inheritanceLoss = 0.0;
        for(int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
            final int[] sampleIndices = trioSampleIndices[trioIndex];
            final int fatherSampleIndex = sampleIndices[0];
            final int motherSampleIndex = sampleIndices[1];
            final int childSampleIndex = sampleIndices[2];
            final int fatherAc = sampleAlleleCounts[fatherSampleIndex];
            final int motherAc = sampleAlleleCounts[motherSampleIndex];
            final int childAc = sampleAlleleCounts[childSampleIndex];
            final Integer fatherPredictionIndex = sampleIndexToPredictionIndex[fatherSampleIndex];
            final Integer motherPredictionIndex = sampleIndexToPredictionIndex[motherSampleIndex];
            final Integer childPredictionIndex = sampleIndexToPredictionIndex[childSampleIndex];
            final float pFather = fatherPredictionIndex == null ? 1F : pSample[fatherPredictionIndex];
            final float pMother = motherPredictionIndex == null ? 1F : pSample[motherPredictionIndex];
            final float pChild = childPredictionIndex == null ? 1F : pSample[childPredictionIndex];

            final double pMendelian = FastMath.max(
                getPMendelian(fatherAc, motherAc, childAc, pFather, pMother, pChild, d1PMendelianDPSample),
                PROB_EPS
            );
            inheritanceLoss -= FastMath.log(pMendelian);
            if(fatherPredictionIndex != null) {
                final double d1 = -d1PMendelianDPSample[0] / pMendelian;
                setLossDerivs(pFather, d1, d1 * d1, inheritanceDerivWeight, d1PSample, d2PSample,
                              fatherSampleIndex);
                if(!(goodSampleIndices.contains(fatherSampleIndex) || badSampleIndices.contains(fatherSampleIndex))) {
                    // add fatherPredictionIndex to stream of indices that have a prior that the upper 50% should be
                    // "good"
                    priorSampleIndices.add(fatherPredictionIndex);
                }
            }
            if(motherPredictionIndex != null) {
                final double d1 = -d1PMendelianDPSample[1] / pMendelian;
                setLossDerivs(pMother, d1, d1 * d1, inheritanceDerivWeight, d1PSample, d2PSample,
                              motherSampleIndex);
                if(!(goodSampleIndices.contains(motherSampleIndex) || badSampleIndices.contains(motherSampleIndex))) {
                    // add motherPredictionIndex to stream of indices that have a prior that the upper 50% should be
                    // "good"
                    priorSampleIndices.add(motherPredictionIndex);
                }
            }
            if(childPredictionIndex != null) {
                final double d1 = -d1PMendelianDPSample[2] / pMendelian;
                setLossDerivs(pChild, d1, d1 * d1, inheritanceDerivWeight, d1PSample, d2PSample,
                              childSampleIndex);
                if(!(goodSampleIndices.contains(childSampleIndex) || badSampleIndices.contains(childSampleIndex))) {
                    // add childPredictionIndex to stream of indices that have a prior that the upper 50% should be
                    // "good"
                    priorSampleIndices.add(childPredictionIndex);
                }

            }
        }
        inheritanceLoss /= numTrios;

        final int numTruth = FastMath.max(1, goodSampleIndices.size() + badSampleIndices.size());
        final double truthDerivWeight = variantTruthWeight * pSample.length / numTruth;
        double truthLoss = 0.0;
        for(final int goodSampleIndex : goodSampleIndices) {
            final Integer predictionIndex = sampleIndexToPredictionIndex[goodSampleIndex];
            if(predictionIndex != null) {
                final double pTruthSample = pSample[predictionIndex];
                final double pGood = FastMath.max(pTruthSample, PROB_EPS);
                truthLoss -= FastMath.log(pGood);
                final double d1 = -1.0 / pGood;
                setLossDerivs(pTruthSample, d1, d1 * d1, truthDerivWeight, d1PSample, d2PSample,
                              goodSampleIndex);
            }
        }
        for(final int badSampleIndex : badSampleIndices) {
            final Integer predictionIndex = sampleIndexToPredictionIndex[badSampleIndex];
            if(predictionIndex != null) {
                final double pTruthSample = pSample[predictionIndex];
                final double pBad = FastMath.max(1.0 - pTruthSample, PROB_EPS);
                truthLoss -= FastMath.log(pBad);
                final double d1 = 1.0 / pBad;
                setLossDerivs(pTruthSample, d1, d1 * d1, truthDerivWeight, d1PSample, d2PSample,
                              badSampleIndex);
            }
        }
        truthLoss /= numTruth;

        return new Loss(inheritanceLoss * variantInheritanceWeight, truthLoss * variantTruthWeight);
    }

    protected Loss getLoss(final float[] pSample, final float[] d1PSample, final float[] d2PSample,
                           final int[] variantIndices) {

        // want average loss per variant, weighted so that each type of propertyBin has equal weight
        final int[] propertyBinCounts = new int[numPropertyBins];
        Arrays.stream(variantIndices)
            .forEach(i -> ++propertyBinCounts[propertyBins[i]]);

        final double[] d1PSampleDouble = new double[numSamples];
        final double[] d2PSampleDouble = new double[numSamples];
        final Integer[] sampleIndexToPredictionIndex = new Integer[numSamples];
        final List<Integer> priorSampleIndices = new ArrayList<>(pSample.length);
        Loss loss = new Loss(0F, 0F);
        int startPredictionIndex = 0;
        for(final int variantIndex : variantIndices) {
            final int propertyBin = propertyBins[variantIndex];
            final double variantInheritanceWeight = propertyBinInheritanceWeights[propertyBin]
                                                  / propertyBinCounts[propertyBin];
            final double variantTruthWeight = propertyBinTruthWeights[propertyBin]
                    / propertyBinCounts[propertyBin];
            loss = Loss.add(
                loss,
                getLoss(pSample, d1PSampleDouble, d2PSampleDouble, startPredictionIndex, sampleIndexToPredictionIndex,
                        variantInheritanceWeight, variantTruthWeight, variantIndex, priorSampleIndices)
            );
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                final Integer predictionIndex = sampleIndexToPredictionIndex[sampleIndex];
                if(predictionIndex != null) {
                    d1PSample[startPredictionIndex] = (float)d1PSampleDouble[sampleIndex];
                    d2PSample[startPredictionIndex] = (float)d2PSampleDouble[sampleIndex];
                    ++startPredictionIndex;
                }
                d1PSampleDouble[sampleIndex] = 0.0;
                d2PSampleDouble[sampleIndex] = 0.0;
            }
        }
        if(startPredictionIndex != pSample.length) {
            throw new GATKException("Didn't actually handle loss of all pSample: final startPredictionIndex = " +
                    startPredictionIndex + ", pSample.length = " + pSample.length);
        }

        // get loss from prior that at least half of the sample-variants should be confidently correct
        //   (loss = -log(p) for the 50% of largest p)
        final int numPenalizeFail = (int)FastMath.round( optimalProportionOfSampleVariantsPassing
                                                            * priorSampleIndices.size());
        final int numAllowFail = priorSampleIndices.size() - numPenalizeFail;
        final Iterator<Integer> bigProbIndices = priorSampleIndices.stream()
                .map(i -> new AbstractMap.SimpleImmutableEntry<>(i, pSample[i]))
                .sorted(Map.Entry.comparingByValue())
                .skip(numAllowFail)
                .mapToInt(AbstractMap.SimpleImmutableEntry::getKey)
                .iterator();
        double priorLoss = 0.0;
        final double priorWeight = 1.0 / numPenalizeFail;
        final double priorDerivWeight = priorWeight * (1.0 - truthWeight) * numPenalizeFail;
        while(bigProbIndices.hasNext()) {
            final int bigProbIndex = bigProbIndices.next();
            final double bigProb = pSample[bigProbIndex];
            priorLoss -= FastMath.log(bigProb);
            final double d1 = -1.0 / bigProb;
            setLossDerivs(bigProb, d1, d1 * d1, priorDerivWeight, d1PSample, d2PSample, bigProbIndex);
        }

        priorLoss *= priorWeight;

        return Loss.add(loss, new Loss(priorLoss, 0F));
    }


    private boolean isMendelianKeepHomvar(final int fatherAc, final int motherAc, final int childAc) {
        // child allele counts should not exhibit de-novo mutations nor be missing inherited homvar
        final int maxAc = (fatherAc > 0 ? 1 : 0) + (motherAc > 0 ? 1 : 0);
        return childAc <= maxAc;
    }


    static protected class FilterSummary {
        final int minGq;
        final long numDiscoverable;
        final long numPassed;
        final long numMendelian;
        final long numTruePositive;
        final long numFalsePositive;
        final long numFalseNegative;
        final String label;

        FilterSummary(final int minGq,
                      final long numDiscoverable, final long numPassed, final long numMendelian,
                      final long numTruePositive, final long numFalsePositive, final long numFalseNegative,
                      final String label) {
            this.minGq = minGq;
            this.numDiscoverable = numDiscoverable;
            this.numPassed = numPassed;
            this.numMendelian = numMendelian;
            this.numTruePositive = numTruePositive;
            this.numFalsePositive = numFalsePositive;
            this.numFalseNegative = numFalseNegative;
            this.label = label;
        }

        static final FilterSummary EMPTY = new FilterSummary(
            Integer.MIN_VALUE, 0L, 0L, 0L, 0L,
            0L, 0L, null
        );

        boolean isEmpty() { return numDiscoverable <= 0; }
        boolean isNotEmpty() { return numDiscoverable > 0; }

        FilterSummary setMinGq(final int newMinGq) {
            return new FilterSummary(newMinGq, numDiscoverable, numPassed, numMendelian,
                                     numTruePositive, numFalsePositive, numFalseNegative, label);
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
                label == null ? other.label : label
            );
        }

        FilterSummary subtract(final FilterSummary other) {
            return new FilterSummary(
                minGq,
                numDiscoverable - other.numDiscoverable, numPassed - other.numPassed,
                numMendelian - other.numMendelian, numTruePositive - other.numTruePositive,
                numFalsePositive - other.numFalsePositive, numFalseNegative - other.numFalseNegative,
                label == null ? other.label : label
            );
        }

        @Override
        public String toString() {
            return label + ":{minGq:" + minGq + ",  numDiscoverable:" + numDiscoverable + ", numPassed:" + numPassed
                   + ", numMendelian:" + numMendelian + ", numTruePositive:" + numTruePositive
                   + ", numFalsePositive:" + numFalsePositive + ", numFalseNegative: " + numFalseNegative + "}";
        }
    }


    protected FilterSummary getFilterSummary(final int minGq, final int variantIndex, final String label) {
        final long numDiscoverable = maxDiscoverableMendelianAc[variantIndex];
        final int[][] variantAlleleCounts = getTrioAlleleCountsMatrix(variantIndex);
        final int[][] variantGenotypeQualities = getTrioGenotypeQualitiesMatrix(variantIndex);
        final int [] goodGqs = goodVariantGqs == null ? null : goodVariantGqs.getOrDefault(variantIndex, null);
        final int [] badGqs = badVariantGqs == null ? null : badVariantGqs.getOrDefault(variantIndex, null);
        return getFilterSummary(minGq, numDiscoverable, variantAlleleCounts, variantGenotypeQualities,
                                goodGqs, badGqs, label);
    }

    protected BinnedFilterSummaries getBinnedFilterSummary(final int minGq, final int variantIndex) {
        final int propertyBin = propertyBins[variantIndex];
        return new BinnedFilterSummaries(getFilterSummary(minGq, variantIndex, propertyBinLabels[propertyBin]),
                                         propertyBin, numPropertyBins);
    }

    protected FilterSummary getFilterSummary(final int minGq, final long numDiscoverable,
                                             final int[][] alleleCounts, final int[][] genotypeQualities,
                                             final int[] goodGqs, final int[] badGqs, final String label) {
        long numPassed = 0;
        long numMendelian = 0;
        long numTruePositives = 0;
        long numFalsePositives = 0;
        long numFalseNegatives = 0;
        if(goodGqs != null) {
            for (final int gq : goodGqs) {
                if(gq >= minGq) {
                    ++numTruePositives;
                } else {
                    ++numFalseNegatives;
                }
            }
        }
        if(badGqs != null) {
            for(final int gq : badGqs) {
                if(gq >= minGq) {
                    ++numFalsePositives;
                }
            }
        }
        if(keepHomvar) {
            for (int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
                final int[] trioAc = alleleCounts[trioIndex];
                final int numFilterableTrio =
                        (trioAc[0] == 1 ? 1 : 0) + (trioAc[1] == 1 ? 1 : 0) + (trioAc[2] == 1 ? 1 : 0);

                if (numFilterableTrio == 0) {
                    continue;
                }
                final int[] trioGq = genotypeQualities[trioIndex];
                final int fatherAc = (trioAc[0] == 1 && trioGq[0] < minGq) ? 0 : trioAc[0];
                final int motherAc = (trioAc[1] == 1 && trioGq[1] < minGq) ? 0 : trioAc[1];
                final int childAc = (trioAc[2] == 1 && trioGq[2] < minGq) ? 0 : trioAc[2];
                // Note that we only consider an allele to have "passed" if it was in principle filterable:
                final int numPassedTrio =
                        (fatherAc == 1 ? 1 : 0) + (motherAc == 1 ? 1 : 0) + (childAc == 1 ? 1 : 0);
                if (numPassedTrio > 0) {
                    numPassed += numPassedTrio;
                    if (isMendelianKeepHomvar(fatherAc, motherAc, childAc)) {
                        numMendelian += numPassedTrio;
                    }
                }
            }
        } else {
            for (int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
                final int[] trioAc = alleleCounts[trioIndex];
                final int numFilterableTrio = trioAc[0] + trioAc[1] + trioAc[2];

                if (numFilterableTrio == 0) {
                    continue;
                }
                final int[] trioGq = genotypeQualities[trioIndex];
                final int fatherAc = trioGq[0] < minGq ? 0 : trioAc[0];
                final int motherAc = trioGq[1] < minGq ? 0 : trioAc[1];
                final int childAc = trioGq[2] < minGq ? 0 : trioAc[2];
                final int numPassedTrio = fatherAc + motherAc + childAc;
                if (numPassedTrio > 0) {
                    numPassed += numPassedTrio;
                    if (isMendelian(fatherAc, motherAc, childAc)) {
                        numMendelian += numPassedTrio;
                    }
                }
            }
        }
        return new FilterSummary(minGq, numDiscoverable, numPassed, numMendelian,
                                 numTruePositives, numFalsePositives, numFalseNegatives, label);
    }

    private int getFilteredAlleleCounts(final int sampleIndex, final int[] sampleAlleleCounts,
                                        final Integer[] sampleIndexToPassIndex, final boolean[] passFilter) {
        final int vcfAlleleCounts = sampleAlleleCounts[sampleIndex];
        final Integer passIndex = sampleIndexToPassIndex[sampleIndex];
        return passIndex == null || passFilter[passIndex] ? vcfAlleleCounts : 0;
    }

    protected FilterSummary getFilterSummary(final boolean[] samplePasses, final int variantIndex, final String label) {
        long numPassed = 0;
        long numMendelian = 0;
        long numTruePositives = 0;
        long numFalsePositives = 0;
        long numFalseNegatives = 0;

        final Set<Integer> goodSampleIndices = goodSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet());
        final Set<Integer> badSampleIndices = badSampleVariantIndices.getOrDefault(variantIndex, Collections.emptySet());
        final int[] sampleAlleleCounts = sampleVariantAlleleCounts.get(variantIndex);
        final Integer[] sampleIndexToPredictionIndex = new Integer[numSamples];
        setSampleIndexToPredictionIndex(sampleIndexToPredictionIndex, sampleAlleleCounts,
                                        goodSampleIndices, badSampleIndices, 0);
        for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
            final Integer predictionIndex = sampleIndexToPredictionIndex[sampleIndex];
            if(predictionIndex == null) {
                continue;
            }
            if(goodSampleIndices.contains(sampleIndex)) {
                if(samplePasses[predictionIndex]) {
                    ++numTruePositives;
                } else {
                    ++numFalseNegatives;
                }
            } else if(badSampleIndices.contains(sampleIndex) && samplePasses[predictionIndex]) {
                ++numFalsePositives;
            }
        }

        for(int trioIndex = 0; trioIndex < numTrios; ++trioIndex) {
            final int[] sampleIndices = trioSampleIndices[trioIndex];
            final int fatherIndex = sampleIndices[0];
            final int motherIndex = sampleIndices[1];
            final int childIndex = sampleIndices[2];
            final boolean fatherIsNotFilterable =  sampleIndexToPredictionIndex[fatherIndex] == null;
            final boolean motherIsNotFilterable =  sampleIndexToPredictionIndex[motherIndex] == null;
            final boolean childIsNotFilterable =  sampleIndexToPredictionIndex[childIndex] == null;
            if(fatherIsNotFilterable && motherIsNotFilterable && childIsNotFilterable) {
                continue;
            }
            final int fatherAc = getFilteredAlleleCounts(fatherIndex, sampleAlleleCounts, sampleIndexToPredictionIndex,
                                                         samplePasses);
            final int motherAc = getFilteredAlleleCounts(motherIndex, sampleAlleleCounts, sampleIndexToPredictionIndex,
                                                         samplePasses);
            final int childAc = getFilteredAlleleCounts(childIndex, sampleAlleleCounts, sampleIndexToPredictionIndex,
                                                        samplePasses);
            // Note that we only consider an allele to have "passed" if it was in principle filterable:
            final int numPassedTrio = (fatherIsNotFilterable ? 0 : fatherAc) +
                    (motherIsNotFilterable ? 0 : motherAc) +
                    (childIsNotFilterable ? 0 : childAc);
            if(numPassedTrio > 0) {
                numPassed += numPassedTrio;
                if(isMendelian(fatherAc, motherAc, childAc)) {
                    numMendelian += numPassedTrio;
                }
            }
        }

        final long numDiscoverable = maxDiscoverableMendelianAc[variantIndex];
        return new FilterSummary(minGQ, numDiscoverable, numPassed, numMendelian,
                                 numTruePositives, numFalsePositives, numFalseNegatives, label);
    }

    protected BinnedFilterSummaries getBinnedFilterSummary(final boolean[] samplePasses, final int variantIndex) {
        final int propertyBin = propertyBins[variantIndex];
        return new BinnedFilterSummaries(getFilterSummary(samplePasses, variantIndex, propertyBinLabels[propertyBin]),
                                         propertyBin, numPropertyBins);
    }


    protected BinnedFilterSummaries getBinnedFilterSummary(final int minGq, final long numDiscoverable,
                                             final int[][] alleleCounts, final int[][] genotypeQualities,
                                             final int[] goodGqs, final int[] badGqs, final int propertyBin) {
        return new BinnedFilterSummaries(
            getFilterSummary(minGq, numDiscoverable, alleleCounts, genotypeQualities, goodGqs, badGqs,
                             propertyBinLabels[propertyBin]),
            propertyBin, numPropertyBins
        );
    }

    static protected double getInheritanceF1(final long numDiscoverableMendelian, final long numMendelian, final long numPassed) {
        // calculate f1 score:
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        //     recall = numMendelian / numDiscoverableMendelian
        //     precision = numMendelian / numPassed
        //     -> f1 = 2.0 * numMendelian / (numNonRef + numPassed)
        // No data to filter? -> no loss -> f1 = 1.0
        final double f1 = (numDiscoverableMendelian > 0) ?
                2.0 * (double)numMendelian / (double)(numDiscoverableMendelian + numPassed) :
                1.0;
        if(f1 < 0 || f1 > 1) {
            throw new GATKException("f1 out of range(f1=" + f1 + "). numDiscoverableMendelian=" + numDiscoverableMendelian
                                    + ", numMendelian=" + numMendelian + ", numPassed=" + numPassed);
        }
        return f1;
    }

    static protected double getTruthF1(final long numTruePositive, final long numFalsePositive, final long numFalseNegative) {
        // NOTE: calls are compared against list of true variants and false variants
        // calculate f1 score:
        //     f1 = 2.0 / (1.0 / recall + 1.0 / precision)
        //     recall = numTruePositive / (numTruePositive + numFalseNegative)
        //     precision = numTruePositive / (numTruePositive + numFalsePositive)
        //     -> f1 = 2.0 * numTruePositive / (2 * numTruePositive + numFalseNegative + numFalsePositive)
        // No data to filter? -> no loss -> f1 = 1.0
        final double f1 = numFalseNegative + numFalsePositive > 0 ?
            2.0 * numTruePositive / (2.0 * numTruePositive + numFalseNegative + numFalsePositive) :
            1.0;
        if(f1 < 0 || f1 > 1) {
            // this should never happen => negative numbers as inputs, i.e. a bug in accumulation of raw stats
            throw new GATKException("f1 out of range(f1=" + f1 + "). numTruePositive=" + numTruePositive
                    + ", numFalsePositive=" + numFalsePositive + ", numFalseNegative=" + numFalseNegative);
        }
        return f1;
    }

    static protected class Loss implements Comparable<Loss> {
        final double inheritanceLoss;
        final double truthLoss;
        final String label;
        final static double truthWeight = MinGqVariantFilterBase.truthWeight;
        final static double inheritanceWeight = 1.0 - truthWeight;

        Loss(final double inheritanceLoss, final double truthLoss, final String label) {
            this.inheritanceLoss = inheritanceLoss;
            this.truthLoss = truthLoss;
            this.label = label;
        }

        Loss(final double inheritanceLoss, final double truthLoss) {
            this.inheritanceLoss = inheritanceLoss;
            this.truthLoss = truthLoss;
            this.label = null;
        }

        Loss(final FilterSummary filterSummary) {
            this(
                    1.0 - getInheritanceF1(filterSummary.numDiscoverable, filterSummary.numMendelian, filterSummary.numPassed),
                    1.0 - getTruthF1(filterSummary.numTruePositive, filterSummary.numFalsePositive, filterSummary.numFalseNegative),
                    filterSummary.label
            );
        }

        Loss(final Loss copyLoss) {
            this(copyLoss.inheritanceLoss, copyLoss.truthLoss, copyLoss.label);
        }

        Loss(final Loss copyLoss, final String label) {
            this(copyLoss.inheritanceLoss, copyLoss.truthLoss, label);
        }

        Loss(final Collection<FilterSummary> filterSummaries) {
            this(
                filterSummaries.stream()
                    .map(Loss::new)
                    .reduce(Loss::add).orElse(EMPTY)
                    .divide(max((int)filterSummaries.stream()
                                        .filter(FilterSummary::isNotEmpty).count(),
                             1) ),
                filterSummaries.stream()
                    .filter(FilterSummary::isNotEmpty)
                    .map(Loss::new)
                    .map(Loss::toString)
                    .collect(Collectors.joining("\n"))
                    + "\noverall: "
            );
        }

        Loss(final float[] pSampleVariantIsGood, final boolean[] sampleVariantIsGood) {
            final double balancedLoss = IntStream.range(0, pSampleVariantIsGood.length)
                .mapToDouble(
                    idx -> -FastMath.log(
                        sampleVariantIsGood[idx] ? (double)pSampleVariantIsGood[idx] : 1.0 - (double)pSampleVariantIsGood[idx]
                    )
                )
                .sum();
            this.inheritanceLoss = balancedLoss;
            this.truthLoss = balancedLoss;
            this.label = null;
        }

        static Loss add(final Loss lossA, final Loss lossB) {
            return new Loss(lossA.inheritanceLoss + lossB.inheritanceLoss,
                            lossA.truthLoss + lossB.truthLoss,
                            lossA.label == null ? lossB.label : lossA.label);
        }

        Loss divide(final int scale) {
            return new Loss(inheritanceLoss / scale, truthLoss / scale, label);
        }

        Loss multiply(final int scale) {
            return new Loss(inheritanceLoss * scale, truthLoss * scale, label);
        }

//        public boolean ge(final Loss other) {
//            return (inheritanceLoss >= other.inheritanceLoss) && (truthLoss >= other.truthLoss);
//        }
//
//        public boolean gt(final Loss other) {
//            return this.ge(other) && (inheritanceLoss > other.inheritanceLoss) || (truthLoss > other.truthLoss);
//        }
//
//        public boolean le(final Loss other) {
//            return (inheritanceLoss <= other.inheritanceLoss) && (truthLoss <= other.truthLoss);
//        }
//
//        public boolean lt(final Loss other) {
//            return this.le(other) && (inheritanceLoss < other.inheritanceLoss) || (truthLoss < other.truthLoss);
//        }

        public boolean ge(final Loss other) {
            return toDouble() >= other.toDouble();
        }

        public boolean gt(final Loss other) {
            return toDouble() > other.toDouble();
        }

        public boolean le(final Loss other) {
            return toDouble() <= other.toDouble();
        }

        public boolean lt(final Loss other) {
            return toDouble() < other.toDouble();
        }

        @Override
        public int compareTo(final @NotNull Loss other) {
            if(this.lt(other)) {
                return -1;
            } else if(this.gt(other)) {
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

        double toDouble() { return inheritanceWeight * inheritanceLoss + truthWeight * truthLoss; }
        float toFloat() { return (float)toDouble(); }

        static final Loss POSITIVE_INFINITY = new Loss(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        static final Loss NaN = new Loss(Double.NaN, Double.NaN);
        static final Loss EMPTY = new Loss(FilterSummary.EMPTY);
    }

    static protected class BinnedFilterSummaries extends ArrayList<FilterSummary> {
        private static final long serialVersionUID = 0;

        BinnedFilterSummaries(final Collection<FilterSummary> filterSummaries) {
            super(filterSummaries);
        }

        BinnedFilterSummaries(final FilterSummary filterSummary, final int binIndex, final int numBins) {
            this(
                IntStream.range(0, numBins)
                    .mapToObj(i -> i == binIndex ? filterSummary : FilterSummary.EMPTY)
                    .collect(Collectors.toList())
            );
        }

        BinnedFilterSummaries add(final BinnedFilterSummaries other) {
            return this.equals(BinnedFilterSummaries.EMPTY) ?
                other :
                other.equals(BinnedFilterSummaries.EMPTY) ?
                    this :
                    new BinnedFilterSummaries(
                        IntStream.range(0, this.size()).mapToObj(
                            i -> this.get(i).add(other.get(i))).collect(Collectors.toList()
                        )
                    );
        }

        BinnedFilterSummaries subtract(final BinnedFilterSummaries other) {
            return other.equals(BinnedFilterSummaries.EMPTY) ?
                this :
                this.equals(BinnedFilterSummaries.EMPTY) ?
                    new BinnedFilterSummaries(
                            other.stream()
                                .map(FilterSummary.EMPTY::subtract)
                                .collect(Collectors.toList())
                    ) :
                    new BinnedFilterSummaries(
                            IntStream.range(0, this.size()).mapToObj(
                                    i -> this.get(i).subtract(other.get(i))).collect(Collectors.toList()
                            )
                    );
        }

        int getMinGq() {
            // Different bins may be empty (set to Integer.MIN_VALUE)
            // Thus minGq is either inconsistent (summed over multiple variants) or the maximum over all bins.
            return this.stream().mapToInt(filterSummary -> filterSummary.minGq).max().orElseThrow(RuntimeException::new);
        }

        BinnedFilterSummaries setMinGq(final int newMinGq) {
            return new BinnedFilterSummaries(
                this.stream().map(filterSummary -> filterSummary.setMinGq(newMinGq)).collect(Collectors.toList())
            );
        }

        static BinnedFilterSummaries EMPTY = new BinnedFilterSummaries(FilterSummary.EMPTY, 0, 1);
    }

    static protected class FilterQuality extends BinnedFilterSummaries implements Comparable<FilterQuality> {
        private static final long serialVersionUID = 0;
        final Loss loss;

        FilterQuality(final Collection<FilterSummary> filterSummaries) {
            super(filterSummaries);
            this.loss = new Loss(filterSummaries);
        }

        FilterQuality(final FilterSummary filterSummary) {
            this(Collections.singletonList(filterSummary));
        }


        @Override
        public int compareTo(final @NotNull FilterQuality other) { return this.loss.compareTo(other.loss); }

        static final FilterQuality EMPTY = new FilterQuality(FilterSummary.EMPTY);
    }

    protected FilterQuality getOptimalVariantMinGq(final int variantIndex, final BinnedFilterSummaries backgroundFilterSummary) {
        final IntStream candidateMinGqs = getCandidateMinGqs(variantIndex);
        if(candidateMinGqs == null) {
            // minGq doesn't matter for this row, so return previous optimal filter or trivial filter
            return backgroundFilterSummary == null ? FilterQuality.EMPTY : new FilterQuality(backgroundFilterSummary);
        }

        final long numDiscoverable = maxDiscoverableMendelianAc[variantIndex];
        final int[][] variantAlleleCounts = getTrioAlleleCountsMatrix(variantIndex);
        final int[][] variantGenotypeQualities = getTrioGenotypeQualitiesMatrix(variantIndex);
        final int [] goodGqs = goodVariantGqs == null ? null : goodVariantGqs.getOrDefault(variantIndex, null);
        final int [] badGqs = badVariantGqs == null ? null : badVariantGqs.getOrDefault(variantIndex, null);
        final int propertyBin = propertyBins[variantIndex];
        if(backgroundFilterSummary == null) {
            // doing optimization only considering loss of each individual variant
            return candidateMinGqs.parallel()
                    .mapToObj(minGq -> new FilterQuality(
                            getBinnedFilterSummary(minGq, numDiscoverable, variantAlleleCounts, variantGenotypeQualities,
                                                   goodGqs, badGqs, propertyBin)
                        ))
                    .min(FilterQuality::compareTo)
                    .orElseThrow(RuntimeException::new);
        } else {
            // doing optimization considering overall loss
            return candidateMinGqs
                    .parallel()
                    .mapToObj(minGq -> new FilterQuality(
                            getBinnedFilterSummary(minGq, numDiscoverable, variantAlleleCounts, variantGenotypeQualities,
                                                   goodGqs, badGqs, propertyBin)
                            .add(backgroundFilterSummary)
                    ))
                    .min(FilterQuality::compareTo)
                    .orElseThrow(RuntimeException::new);
        }
    }

    protected BinnedFilterSummaries getBackgroundFilterSummary(
            final int[] optimizingIndices, final int[] trainingIndices,
            int[] minGqs
    ) {
        // Make empty summary
        BinnedFilterSummaries filterSummary = BinnedFilterSummaries.EMPTY;
        if(trainingIndices == null || trainingIndices.length == 0 || trainingIndices == optimizingIndices) {
            // Not using any background, return trivial background summary
            return filterSummary;
        }
        if(minGqs == null) {
            throw new GATKException("If using non-trivial trainingIndices, must pass non-null minGqs with length equal to trainingIndices.length");
        }

        int optimizeIndex = 0;
        int nextRow = optimizingIndices == null ? Integer.MAX_VALUE : optimizingIndices[optimizeIndex];
        if(nextRow < trainingIndices[0]) {
            throw new GATKException("optimizingIndices start before training set");
        }
        // Add results at training indices not in the eval set
        for (int i = 0; i < trainingIndices.length; ++i) {
            final int trainingIndex = trainingIndices[i];
            if (nextRow > trainingIndex) {
                // This index is not in the eval set, add it to background
                filterSummary = filterSummary.add(getBinnedFilterSummary(minGqs[i], trainingIndex));
            } else {
                // This index is in the optimize set, skip it and point to next optimize index
                ++optimizeIndex;
                nextRow = optimizeIndex < optimizingIndices.length ?
                    optimizingIndices[optimizeIndex] :
                    Integer.MAX_VALUE;
            }
        }

        return filterSummary;
    }

    /**
     * Optimize minGq as a constant value on optimizeIndices,
     * @param optimizeIndices
     * @param trainingIndices
     * @param minGqs
     * @return
     */
    protected FilterQuality getOptimalMinGq(final int[] optimizeIndices,
                                            final int[] trainingIndices, final int[] minGqs) {
        if(optimizeIndices.length == 0) {
            throw new GATKException("Can't get optimalMinGq from empty rowIndices");
        }
        IntStream candidateMinGqs = getCandidateMinGqs(optimizeIndices);
        if(candidateMinGqs == null) {
            // minGq doesn't matter for these rows, just return something
            candidateMinGqs = IntStream.of(Integer.MIN_VALUE);
        }

        final BinnedFilterSummaries backgroundFilterSummary = getBackgroundFilterSummary(
                optimizeIndices, trainingIndices, minGqs
        );

        // For each candidate minGq:
        //      add backgroundFilterSummary to FilterSummaries from optimizeIndices
        //      convert to FilterQuality (calculate Loss)
        // Return FilterQuality corresponding to best Loss
        return candidateMinGqs
                .parallel()
                .mapToObj(
                    candidateMinGq -> new FilterQuality(
                        Arrays.stream(optimizeIndices)
                            .mapToObj(evalIndex -> getBinnedFilterSummary(candidateMinGq, evalIndex))
                            .reduce(backgroundFilterSummary.setMinGq(candidateMinGq), BinnedFilterSummaries::add)
                    )
                )
                .min(FilterQuality::compareTo)
                .orElseThrow(() -> new GATKException("Could not find optimal minGq filter. This is a bug."));
    }

    private void setMaxDiscoverableMendelianAc() {
        maxDiscoverableMendelianAc = new int [numVariants];
        numDiscoverableMendelianAc = 0;
        for(int variantIndex = 0; variantIndex < numVariants; ++variantIndex) {
            final int[][] variantAlleleCounts = getTrioAlleleCountsMatrix(variantIndex);
            final int[][] variantGenotypeQualities = getTrioGenotypeQualitiesMatrix(variantIndex);
            final IntStream candidateMinGq = getCandidateMinGqs(variantAlleleCounts, variantGenotypeQualities, null);

            maxDiscoverableMendelianAc[variantIndex] = candidateMinGq == null ? 0 :
                candidateMinGq.parallel().map(
                        minGq -> (int)getFilterSummary(minGq, 0, variantAlleleCounts, variantGenotypeQualities, null, null, null).numMendelian
                ).max().orElse(0);
            numDiscoverableMendelianAc += maxDiscoverableMendelianAc[variantIndex];
        }
    }

    private void setOptimalProportionOfSampleVariantsPassing() {
        int numPass = 0;
        int numFilterable = 0;

        for(int variantIndex = 0; variantIndex < numVariants; ++variantIndex) {
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(getSampleVariantIsFilterable(variantIndex, sampleIndex)) {
                    ++numFilterable;
                    final int gq = sampleVariantGenotypeQualities.get(variantIndex)[sampleIndex];
                    if(gq >= perVariantOptimalMinGq[variantIndex]) {
                        ++numPass;
                    }
                }
            }
        }
        optimalProportionOfSampleVariantsPassing = numPass / (double)numFilterable;
        if(printProgress > 0) {
            System.out.format("Optimal proportion of variants passing: %.3f\n", optimalProportionOfSampleVariantsPassing);
        }
    }

    private void setPerVariantOptimalMinGq() {
        // Get initial optimal filter qualities, optimizing each variant separately
        // Collect total summary stats, store min GQ
        perVariantOptimalMinGq = new int[numVariants];
        BinnedFilterSummaries overallSummary = BinnedFilterSummaries.EMPTY;
        final BinnedFilterSummaries[] filterSummaries = new BinnedFilterSummaries[numVariants];
        for(int variantIndex = 0; variantIndex < numVariants; ++variantIndex) {
            final FilterQuality greedyFilter = getOptimalVariantMinGq(variantIndex, null);
            filterSummaries[variantIndex] = greedyFilter;
            overallSummary = overallSummary.add(greedyFilter);
            perVariantOptimalMinGq[variantIndex] = greedyFilter.getMinGq();
        }

        // Iteratively improve filters, optimizing for OVERALL loss
        boolean anyImproved = true;
        int numIterations = 0;
        Loss previousLoss = new Loss(overallSummary);
        while(anyImproved) {
            ++numIterations;
            anyImproved = false;
            for(int variantIndex = 0; variantIndex < numVariants; ++variantIndex) {
                final BinnedFilterSummaries previousFilter = filterSummaries[variantIndex];
                final BinnedFilterSummaries backgroundFilter = overallSummary.subtract(previousFilter);
                final FilterQuality greedyFilter = getOptimalVariantMinGq(variantIndex, backgroundFilter);
                perVariantOptimalMinGq[variantIndex] = greedyFilter.getMinGq();
                filterSummaries[variantIndex] = greedyFilter.subtract(backgroundFilter);
                overallSummary = greedyFilter;

                if(greedyFilter.loss.lt(previousLoss)) {
                    anyImproved = true;
                } else if(greedyFilter.loss.gt(previousLoss)) {
                    throw new GATKException("Loss increased. This is a bug. previousLoss=" + previousLoss
                                            + ", greedyFilter.loss=" + greedyFilter.loss);
                }
                previousLoss = greedyFilter.loss;
            }
            System.out.println("Iteration " + numIterations + ": loss=" + previousLoss);
        }
        displayHistogram("Optimal minGq histogram", Arrays.stream(perVariantOptimalMinGq), true);

        // add up the variant BinnedFilterSummaries,
        // compute the loss for each bin
        // and set weight for truth and inheritance for each bin based on the loss / performance of the overall
        //      optimal filterSummaries
        overallSummary = Arrays.stream(filterSummaries).reduce(BinnedFilterSummaries::add).orElseThrow(RuntimeException::new);
        final List<Loss> optimalBinnedLosses = overallSummary.stream().map(Loss::new).collect(Collectors.toList());
        propertyBinInheritanceWeights = optimalBinnedLosses.stream()
            .mapToDouble(loss -> loss.inheritanceLoss < 0.5 ? 1.0 - 2 * loss.inheritanceLoss : 0.0)
            .toArray();
        propertyBinTruthWeights = optimalBinnedLosses.stream()
                .mapToDouble(loss -> loss.truthLoss < 0.5 ? 1.0 - 2 * loss.truthLoss : 0.0)
                .toArray();
        final double totalInherit = Arrays.stream(propertyBinInheritanceWeights).sum();
        propertyBinInheritanceWeights = Arrays.stream(propertyBinInheritanceWeights).map(w -> w / totalInherit).toArray();
        final double totalTruth = Arrays.stream(propertyBinTruthWeights).sum();
        propertyBinTruthWeights = Arrays.stream(propertyBinTruthWeights).map(w -> w / totalTruth).toArray();


        // compute the proportion of filterable variants that passed the optimal minGQ filter
        final long numFilterable = IntStream.range(0, numVariants)
                .mapToLong(
                        variantIndex -> IntStream.range(0, numSamples)
                                .filter(sampleIndex -> getSampleVariantIsFilterable(variantIndex, sampleIndex))
                                .count()
                )
                .sum();
        final FilterSummary unBinnedSummary = overallSummary.stream().reduce(FilterSummary::add).orElseThrow(RuntimeException::new);
        optimalProportionOfSampleVariantsPassing = unBinnedSummary.numPassed / (double)numFilterable;
        if(printProgress > 0) {
            System.out.format("Optimal proportion of variants passing: %.3f\n", optimalProportionOfSampleVariantsPassing);
        }
    }



    static int[] take(final int[] values, final int[] indices) {
        return Arrays.stream(indices).map(i -> values[i]).toArray();
    }

    static double[] take(final double[] values, final int[] indices) {
        return Arrays.stream(indices).mapToDouble(i -> values[i]).toArray();
    }

    protected int[] getPerVariantOptimalMinGq(final int[] variantIndices) {
        return take(perVariantOptimalMinGq, variantIndices);
    }

    protected long getNumTrainableSamples(final int variantIndex) {
        return IntStream.range(0, getNumSamples())
                .filter(sampleIndex -> getSampleVariantIsTrainable(variantIndex, sampleIndex))
                .count();
    }
    /**
     * Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
     */
    protected long getNumTrainableSampleVariants(final int[] variantIndices) {
        return Arrays.stream(variantIndices)
            .mapToLong(this::getNumTrainableSamples)
            .sum();
    }

    /**
     *
     * @param variantIndices
     * @return
     */
    protected boolean[] getSampleVariantTruth(final int[] variantIndices) {
        final int numRows = (int) getNumTrainableSampleVariants(variantIndices);
        final boolean[] sampleVariantTruth = new boolean[numRows];

        int flatIndex = 0;
        final PropertiesTable.Property gqMat = propertiesTable.get(VCFConstants.GENOTYPE_QUALITY_KEY);
        for(final int variantIndex : variantIndices) {
            final int minGq = perVariantOptimalMinGq[variantIndex];
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                final int gq = gqMat.getAsInt(variantIndex, sampleIndex);
                sampleVariantTruth[flatIndex] = gq >= minGq;
                ++flatIndex;
            }
        }
        return sampleVariantTruth;
    }

    @SuppressWarnings("SameParameterValue")
    protected void displayHistogram(final String description, final IntStream intStream, boolean binValues) {
        final Map<Integer, Integer> rawValuesMap = new HashMap<>();
        intStream.forEach(gq -> {
            if (rawValuesMap.containsKey(gq)) {
                rawValuesMap.put(gq, 1 + rawValuesMap.get(gq));
            } else {
                rawValuesMap.put(gq, 1);
            }
        });
        if(rawValuesMap.size() == 0) {
            System.out.println(description + ": no data");
        }
        final int minGqValue = rawValuesMap.keySet().stream().min(Integer::compareTo).orElseThrow(RuntimeException::new);
        final int maxGqValue = rawValuesMap.keySet().stream().max(Integer::compareTo).orElseThrow(RuntimeException::new);

        final Map<Integer, Integer> displayValuesMap;
        if(binValues) {
            displayValuesMap = new HashMap<>();
            rawValuesMap.forEach((gq, numGq) -> {
                final int binGq;
                if (gq == 0) {
                    binGq = gq;
                } else {
                    final int magnitude = (int) Math.pow(10.0, Math.floor(Math.log10(Math.abs(gq))));
                    binGq = magnitude * (gq / magnitude);
                }

                if (displayValuesMap.containsKey(binGq)) {
                    displayValuesMap.put(binGq, numGq + displayValuesMap.get(binGq));
                } else {
                    displayValuesMap.put(binGq, numGq);
                }
            });
        } else {
            displayValuesMap = rawValuesMap;
        }
        System.out.println(description + ":");
        System.out.println("min=" + minGqValue + ", max=" + maxGqValue);
        displayValuesMap.keySet()
            .stream()
            .sorted()
            .forEach(minGq -> System.out.println(minGq + ": " + displayValuesMap.get(minGq)));
    }

    @SuppressWarnings("SameParameterValue")
    protected void displayHistogram(final String description, final DoubleStream doubleStream, final int numBins,
                                    final double minValue, final double maxValue) {
        final double safety = 100;
        final double valueEps = safety * FastMath.max(FastMath.ulp(FastMath.abs(minValue)),
                                                      FastMath.ulp(FastMath.abs(maxValue)));
        final double numBinsEps = safety * FastMath.ulp((double)numBins);
        final double valueOffset = minValue - valueEps;
        final double binScale = (numBins - numBinsEps) / (maxValue - valueOffset);
        final long[] valueBins = new long [numBins];
        double[] valueRange = new double[] {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        doubleStream.forEach(value -> {
            final int bin = (int)FastMath.floor((value - valueOffset) * binScale);
            try {
                ++valueBins[bin];
            } catch(ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
                final String errStr = String.format("bad bin: min: %f, max: %f, scale: %f, offset: %f, value: %f\n",
                                                    minValue, maxValue, binScale, valueOffset, value);
                throw new GATKException(errStr, arrayIndexOutOfBoundsException);
            }
            if(value < valueRange[0]) {
                valueRange[0] = value;
            }
            if(value > valueRange[1]) {
                valueRange[1] = value;
            }
        });
        final long numValues = Arrays.stream(valueBins).sum();
        if(valueRange[0] > valueRange[1]) {
            System.out.println(description + ": no data");
        } else {
            System.out.format("%s: %d values\n", description, numValues);
            System.out.format("\tactual range: [%f, %f]\n", valueRange[0], valueRange[1]);
            System.out.println("\t low - high    %");
            double high = minValue;
            for(int bin = 0; bin < numBins; ++bin) {
                double low = high;
                high = low + 1.0 / binScale;
                System.out.format("\t%.2f - %.2f   %.1f\n", low, high, valueBins[bin] * 100.0 / numValues);
            }
        }
    }

    void printDebugInfo() {
        System.out.println("########################################");
        System.out.println("numVariants: " + numVariants);
        System.out.println("numTrios: " + numTrios);
        System.out.println("numProperties: " + getNumProperties());
        System.out.println("index\tpropertyName\tpropertyBaseline\tpropertyScale");
        int idx = 0;
        for(final PropertiesTable.Property property : propertiesTable) {
            System.out.println(idx + "\t" + property.name + "\t" + property.getBaseline() + "\t" + property.getScale());
            ++idx;
        }

        final Map<String, List<String>> labelsEncoding = propertiesTable.getLabelsEncoding();
        for(final Map.Entry<String, List<String>> entry : labelsEncoding.entrySet()) {
            System.out.println(entry.getKey() + ":");
            idx = 0;
            for(final String label : entry.getValue()) {
                System.out.println(idx + "\t" + label);
            }
        }

        final IntStream acStream = sampleVariantAlleleCounts.stream().flatMapToInt(Arrays::stream);
        final IntStream gqStream = sampleVariantGenotypeQualities.stream().flatMapToInt(Arrays::stream);
        displayHistogram("Filterable alleles Gq histogram:",
                            streamFilterableGq(acStream, gqStream),true);

        System.out.println("########################################");
    }

    protected Loss getLoss2(final float[] pSampleVariantGood, final int[] variantIndices) {
        final boolean[] sampleVariantIsGood = getSampleVariantTruth(variantIndices);
        return new Loss(pSampleVariantGood, sampleVariantIsGood);
    }

    protected Loss getLoss(final float[] pSampleVariantGood, final int[] variantIndices) {
        final boolean[][] sampleVariantPasses = new boolean[variantIndices.length][];
        int flatIndex = 0;
        for(int row = 0; row < variantIndices.length; ++row) {
            final int variantIndex = variantIndices[row];
            final boolean[] samplePasses = new boolean[(int)getNumTrainableSamples(variantIndex)];
            sampleVariantPasses[row] = samplePasses;
            for(int col = 0; col < samplePasses.length; ++col) {
                samplePasses[col] = pSampleVariantGood[flatIndex] >= 0.5F;
                ++flatIndex;
            }
        }

        return new Loss(
            IntStream.range(0, variantIndices.length)
                .mapToObj(i -> getBinnedFilterSummary(sampleVariantPasses[i], variantIndices[i]))
                .reduce(BinnedFilterSummaries::add).orElse(BinnedFilterSummaries.EMPTY)
        );
    }

    protected Loss getLoss(final int[] minGq, final int[] variantIndices) {
        if(minGq.length != variantIndices.length) {
            throw new GATKException(
                    "Length of minGq (" + minGq.length + ") does not match length of variantIndices (" + variantIndices.length + ")"
            );
        }
        return new Loss(
            IntStream.range(0, minGq.length)
                .mapToObj(i -> getBinnedFilterSummary(minGq[i], variantIndices[i]))
                .reduce(BinnedFilterSummaries::add).orElse(BinnedFilterSummaries.EMPTY)
        );
    }

    private void setTrainingAndValidationIndices() {
        final int numValidationIndices = (int)round(validationProportion * numVariants);
        final List<Integer> shuffleIndices = IntStream.range(0, numVariants).boxed().collect(Collectors.toList());
        Collections.shuffle(shuffleIndices, randomGenerator);

        validationIndices = shuffleIndices.subList(0, numValidationIndices).stream()
            .sorted().mapToInt(Integer::intValue).toArray();
        trainingIndices = shuffleIndices.subList(numValidationIndices, numVariants).stream()
            .sorted().mapToInt(Integer::intValue).toArray();
    }

    private void saveTrainedModel() {
        try (final OutputStream outputStream = modelFile.getOutputStream()) {
            final OutputStream unclosableOutputStream = new FilterOutputStream(outputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.saveDataEncoding(unclosableOutputStream);
            unclosableOutputStream.write("\n".getBytes());
            saveModel(unclosableOutputStream);
        } catch(IOException ioException) {
            throw new GATKException("Error saving modelFile " + modelFile, ioException);
        }
    }

    private void savePropertiesTable() {
        if(propertiesTableFile == null) {
            return;
        }
        try (final OutputStream outputStream = propertiesTableFile.getOutputStream()) {
            final OutputStream unclosableOutputStream = new FilterOutputStream(outputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.save(unclosableOutputStream);
        } catch(IOException ioException) {
            throw new GATKException("Error saving propertiesTable " + propertiesTableFile, ioException);
        }
    }

    private void loadTrainedModel() {
        if(modelFile == null || !Files.exists(modelFile.toPath())) {
            if(runMode == RunMode.FILTER) {
                throw new UserException("mode=FILTER, but trained model file was not provided.");
            }
            return;
        }
        try (final InputStream inputStream = modelFile.getInputStream()) {
            final InputStream unclosableInputStream = new FilterInputStream(inputStream) {
                @Override
                public void close() {
                    // don't close the stream in one of the subroutines
                }
            };
            propertiesTable.loadDataEncoding(unclosableInputStream);
            loadModel(unclosableInputStream );
        } catch (Exception exception) {
            throw new GATKException("Error loading modelFile " + modelFile + " (malformed file?)", exception);
        }
    }

    protected double getDoubleFromJSON(final Object jsonObject) {
        if(jsonObject instanceof Double) {
            return (Double) jsonObject;
        } else if(jsonObject instanceof BigDecimal) {
            return ((BigDecimal)jsonObject).doubleValue();
        } else {
            throw new GATKException("Unknown conversion to double for " + jsonObject.getClass().getName());
        }
    }

    protected List<String> getStringListFromJSON(final Object jsonObject) {
        return ((JSONArray)jsonObject).stream().map(o -> (String)o).collect(Collectors.toList());
    }

    private byte[] modelCheckpoint = null;

    protected void saveModelCheckpoint() {
        final ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
        saveModel(outputStream);
        modelCheckpoint = outputStream.toByteArray();
    }

    protected void loadModelCheckpoint() {
        final ByteArrayInputStream inputStream = new ByteArrayInputStream(modelCheckpoint);
        loadModel(inputStream);
    }

    protected abstract boolean needsNormalizedProperties();
    protected abstract float predict(final float[] sampleVariantProperties);
    protected abstract void trainFilter();
    protected abstract void saveModel(final OutputStream outputStream);
    protected abstract void loadModel(final InputStream inputStream);

    protected int phredScale(final float p) {
        return (int)FastMath.round(-10 * FastMath.log10(p));
    }
    protected int phredScale(final double p) {
        return (int)FastMath.round(-10 * FastMath.log10(p));
    }
    protected int adjustedGq(final float[] sampleVariantProperties) {
        return phredScale(predict(sampleVariantProperties));
    }

    @Override
    public Object onTraversalSuccess() {
        if(runMode == RunMode.TRAIN) {
            if(numVariants == 0) {
                throw new GATKException("No variants contained in vcf: " + drivingVariantFile);
            }

            if(goodSampleVariants != null) {
                // no longer needed, clear them out
                goodSampleVariants.clear();
                goodSampleVariants = null;
                badSampleVariants.clear();
                badSampleVariants = null;
            }

            setPropertyBins();
            propertiesTable.validateAndFinalize();

            savePropertiesTable();
            setTrainingAndValidationIndices();
            setMaxDiscoverableMendelianAc();
            setPerVariantOptimalMinGq();

            printDebugInfo();

            trainFilter();
            saveTrainedModel();
        } else {
            System.out.println("Filtered " + numFilteredGenotypes + " genotypes in " + numVariants + " variants.");
        }
        return null;
    }
}
