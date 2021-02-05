package org.broadinstitute.hellbender.tools.walkers.sv;

import net.minidev.json.JSONArray;
import net.minidev.json.JSONObject;
import net.minidev.json.JSONValue;
import net.minidev.json.parser.ParseException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.util.*;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.IntStream;

/**
 * Class to manage table with mixed columnar and matrix properties; with different primitive types. Supported types:
 *     boolean, int, long, float, double
 * Conceptually the table is organized as numRows x numProperties x numColumns
 * Properties are stored as Object that can be cast to primitive arrays or array-of-arrays (matrix) in their original
 *     primitive type
 * Rows are extracted as float[] (suitable for machine learning), with one of two options
 *    raw: translated to float type but otherwise unchanged
 *    normalized: shifted so the median is 0, and scaled so the standard deviation (over the central half of the data)
 *                is 1.0. This is done on the fly when extracting rows.
 *                The baseline and scale can be provided when adding a column (to provide consistency between training
 *                and inference data) or computed automatically.
 * This class only requires a lot of code because it's not possible to use generics on primitive types...
 */
class PropertiesTable {
    private static final String PROPERTY_NAMES_KEY = "propertyNames";
    private static final String PROPERTY_BASELINE_KEY = "propertyBaseline";
    private static final String PROPERTY_SCALE_KEY = "propertyScale";

    private static final String BOOLEAN_ARR = "boolean[]";
    private static final String BOOLEAN_MAT = "boolean[][]";
    private static final String INT_ARR = "int[]";
    private static final String INT_MAT = "int[][]";
    private static final String LONG_ARR = "long[]";
    private static final String LONG_MAT = "long[][]";
    private static final String FLOAT_ARR = "float[]";
    private static final String FLOAT_MAT = "float[][]";
    private static final String DOUBLE_ARR = "double[]";
    private static final String DOUBLE_MAT = "double[][]";

    private final List<Object> properties = new ArrayList<>();
    private final List<String> propertyNames = new ArrayList<>();
    private final List<String> propertyClasses = new ArrayList<>();
    private final List<Float> baselines = new ArrayList<>();
    private final List<Float> scales = new ArrayList<>();

    PropertiesTable() {}

    public void addValues(final String propertyName, final boolean[] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final long numTrue = getNumTrue(values);
            final long numFalse = values.length - numTrue;
            baseline = numTrue / (double) values.length;
            scale = numTrue == 0 || numFalse == 0 ?
                    1.0 :
                    FastMath.sqrt(baseline * numFalse / ((double) values.length));
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    protected static long getNumTrue(final boolean[] values) {
        long numTrue = 0;
        for(final boolean value : values) {
            if(value) {
                ++numTrue;
            }
        }
        return numTrue;
    }

    public void addValues(final String propertyName, final boolean[][] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final long numTrue = Arrays.stream(values).mapToLong(this::getNumTrue).sum();
            final long numFalse = values.length - numTrue;
            baseline = numTrue / (double) values.length;
            scale = numTrue == 0 || numFalse == 0 ?
                    1.0 :
                    FastMath.sqrt(baseline * numFalse / ((double) values.length));
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final int[] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).sorted().mapToDouble(x -> x).toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final int[][] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).flatMapToInt(Arrays::stream).sorted()
                .mapToDouble(x -> x).toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final long[] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).sorted().mapToDouble(x -> x).toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final long[][] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).flatMapToLong(Arrays::stream).sorted()
                .mapToDouble(x -> x).toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final float[] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = IntStream.range(0, values.length)
                .mapToDouble(i -> (double)values[i])
                .sorted()
                .toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final float[][] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values)
                .flatMapToDouble(arr -> IntStream.range(0, arr.length).mapToDouble(i -> arr[i]))
                .sorted()
                .toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final double[] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).sorted().toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    public void addValues(final String propertyName, final double[][] values, Double baseline, Double scale) {
        if(baseline == null || scale == null) {
            if(!(baseline == null && scale == null)) {
                throw new IllegalArgumentException("baseline and scale must both be provided or both be omitted");
            }
            final double[] orderedValues = Arrays.stream(values).flatMapToDouble(Arrays::stream).sorted().toArray();
            baseline = getBaselineOrdered(orderedValues);
            scale = getScaleOrdered(orderedValues, baseline);
        }
        addValues(propertyName, (Object)values, baseline, scale);
    }

    protected void addValues(final String propertyName, final Object values, final double baseline, final double scale) {
        // Maintain propertyNames in sorted order, and keep other entries in same order
        final int keySearchIndex = Collections.binarySearch(propertyNames, propertyName);
        if(keySearchIndex >= 0) {
            throw new IllegalArgumentException("PropertiesTable already contains propertyName \"" + propertyName + "\"");
        } else {
            final int insertIndex = -1 - keySearchIndex;
            propertyNames.add(insertIndex, propertyName);
            properties.add(insertIndex, values);
            baselines.add(insertIndex, (float)baseline);
            scales.add(insertIndex, (float)scale);
        }
    }


    protected float getAsFloat(final int keyIndex, final int rowIndex, final int hyperIndex, final boolean normalize) {
        final float rawValue;
        final Object valuesObject = properties.get(keyIndex);
        final String propertyClassName = propertyClasses.get(keyIndex);
        switch(propertyClassName) {
            case BOOLEAN_ARR:
                rawValue = ((boolean[]) valuesObject)[rowIndex] ? 1F : 0F;
                break;
            case INT_ARR:
                rawValue = ((int[])valuesObject)[rowIndex];
                break;
            case LONG_ARR:
                rawValue = ((long[])valuesObject)[rowIndex];
                break;
            case FLOAT_ARR:
                rawValue = ((float[])valuesObject)[rowIndex];
                break;
            case DOUBLE_ARR:
                rawValue = (float)((double[])valuesObject)[rowIndex];
                break;
            case BOOLEAN_MAT:
                rawValue = ((boolean[][]) valuesObject)[rowIndex][hyperIndex] ? 1F : 0F;
                break;
            case INT_MAT:
                rawValue = ((int[][])valuesObject)[rowIndex][hyperIndex];
                break;
            case LONG_MAT:
                rawValue = ((long[][])valuesObject)[rowIndex][hyperIndex];
                break;
            case FLOAT_MAT:
                rawValue = ((float[][])valuesObject)[rowIndex][hyperIndex];
                break;
            case DOUBLE_MAT:
                rawValue = (float)((double[][])valuesObject)[rowIndex][hyperIndex];
                break;
            default:
                throw new IllegalArgumentException(
                    "Error getting " + propertyNames.get(keyIndex) + ": PropertiesTable does not handle values of type \""
                    + propertyClassName + "\""
                );
        }
        return normalize ?
                (rawValue - baselines.get(keyIndex)) / scales.get(keyIndex) :
                rawValue;
    }

    public int fillProperties(final int rowIndex, final float[] outArray, int outIndex, final int hyperIndex,
                              final boolean normalize) {
        for(int keyIndex = 0; keyIndex < propertyNames.size(); ++keyIndex) {
            outArray[outIndex] = getAsFloat(keyIndex, rowIndex, hyperIndex, normalize);
            ++outIndex;
        }
        return outIndex;
    }

    static protected double getBaselineOrdered(final double[] orderedValues) {
        // get baseline as median of values
        return orderedValues.length == 0 ?
                0 :
                orderedValues.length % 2 == 1 ?
                        orderedValues[orderedValues.length / 2] :
                        (orderedValues[orderedValues.length / 2 - 1] + orderedValues[orderedValues.length / 2]) / 2.0;
    }

    static protected double getScaleOrdered(final double[] orderedValues, final double baseline) {
        // get scale as root-mean-square difference from baseline, over central half of data (to exclude outliers)
        switch(orderedValues.length) {
            case 0:
            case 1:
                return 1.0;
            default:
                final int start = orderedValues.length / 4;
                final int stop = 3 * orderedValues.length / 4;
                double scale = 0.0;
                for(int idx = start; idx < stop; ++idx) {
                    scale += (orderedValues[idx] - baseline) * (orderedValues[idx] - baseline);
                }
                return FastMath.max(FastMath.sqrt(scale / (1 + stop - start)), 1.0e-6);
        }
    }

    static abstract class Property {
        protected Float baseline = null;
        protected Float scale = null;
        protected int numRows;

        static final int initialAllocatedRows = 100;
        static final int allocationGrowthScale = 2;

        abstract public float getAsFloat(final int rowIndex, final int hyperIndex);
        abstract protected int getAllocatedRows();
        abstract public void setAllocatedRowsUnguarded(final int numRows);
        abstract protected void assignNextValue(final Object value);
        abstract protected double[] getValuesAsOrderedDoubles();

        public int getNumRows() { return numRows; }

        public void setAllocatedRows(final int numRows) {
            if(numRows <= getNumRows()) {
                if(numRows < getNumRows()) {
                    throw new IllegalArgumentException("Can't shrink allocated memory for property to less than used space.");
                }
            } else {
                setAllocatedRowsUnguarded(numRows);
            }
        }

        public float getAsFloat(final int rowIndex, final int hyperIndex, final boolean normalize) {
            final float rawValue = getAsFloat(rowIndex, hyperIndex);
            return normalize ?
                    (rawValue - baseline) / scale :
                    rawValue;
        }

        public void append(final Object value) {
            if(getNumRows() >= getAllocatedRows()) {
                setAllocatedRows(getNumRows() * allocationGrowthScale);
            }
            assignNextValue(value);
            ++this.numRows;
        }

        protected void calculateScaleAndBaseline() {
            final double[] orderedValues = getValuesAsOrderedDoubles();
            final double baseline = getBaselineOrdered(orderedValues);
            this.baseline = (float)baseline;
            scale = (float)getScaleOrdered(orderedValues, baseline);
        }
    }

    static public class BooleanArrProperty extends Property {
        private boolean[] values = new boolean[initialAllocatedRows];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex] ? 1F : 0F;
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (boolean)value; }
        @Override protected void calculateScaleAndBaseline() {
            // special case for booleans
            final long numTrue = getNumTrue(values);
            final long numFalse = values.length - numTrue;
            final double baseline = numTrue / (double) values.length;
            this.baseline = (float)baseline;
            scale = numTrue == 0 || numFalse == 0 ?
                    1F :
                    (float)FastMath.sqrt(baseline * numFalse / ((double) values.length));
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            // special case for booleans
            throw new GATKException("Method not needed for booleans");
        }
    }

    static public class BooleanMatProperty extends Property {
        private boolean[][] values = new boolean[initialAllocatedRows][];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex] ? 1F : 0F;
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (boolean[])value; }
        @Override protected void calculateScaleAndBaseline() {
            // special case for booleans
            final long numTrue = Arrays.stream(values).mapToLong(PropertiesTable::getNumTrue).sum();
            final long numFalse = values.length - numTrue;
            final double baseline = numTrue / (double) values.length;
            this.baseline = (float)baseline;
            scale = numTrue == 0 || numFalse == 0 ?
                    1F :
                    (float)FastMath.sqrt(baseline * numFalse / ((double) values.length));
        }
        @Override protected double[] getValuesAsOrderedDoubles() {
            // special case for booleans
            throw new GATKException("Method not needed for booleans");
        }
    }

    static public class IntArrProperty extends Property {
        private int[] values = new int[initialAllocatedRows];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (int)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class IntMatProperty extends Property {
        private int[][] values = new int[initialAllocatedRows][];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (int[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).flatMapToInt(Arrays::stream).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class LongArrProperty extends Property {
        private long[] values = new long[initialAllocatedRows];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (long)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class LongMatProperty extends Property {
        private long[][] values = new long[initialAllocatedRows][];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (long[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).flatMapToLong(Arrays::stream).sorted().mapToDouble(x -> x).toArray();
        }
    }

    static public class FloatArrProperty extends Property {
        private float[] values = new float[initialAllocatedRows];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (float)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return IntStream.range(0, values.length).mapToDouble(i -> (double)values[i]).sorted().toArray();
        }
    }

    static public class FloatMatProperty extends Property {
        private float[][] values = new float[initialAllocatedRows][];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return values[rowIndex][hyperIndex];
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (float[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values)
                .flatMapToDouble(arr -> IntStream.range(0, arr.length).mapToDouble(i -> arr[i]))
                .sorted()
                .toArray();
        }
    }

    static public class DoubleArrProperty extends Property {
        private double[] values = new double[initialAllocatedRows];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) { return (float)values[rowIndex]; }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (float)value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).sorted().toArray();
        }
    }

    static public class DoubleMatProperty extends Property {
        private double[][] values = new double[initialAllocatedRows][];

        @Override public float getAsFloat(final int rowIndex, final int hyperIndex) {
            return (float)values[rowIndex][hyperIndex];
        }
        @Override public int getAllocatedRows() { return values.length; }
        @Override public void setAllocatedRowsUnguarded(final int numRows) {
            this.values = Arrays.copyOf(this.values, numRows);
        }
        @Override protected void assignNextValue(final Object value) { this.values[numRows] = (double[])value; }
        @Override protected double[] getValuesAsOrderedDoubles() {
            return Arrays.stream(values).flatMapToDouble(Arrays::stream).sorted().toArray();
        }
    }

    public void saveNormalizationStats(final OutputStream outputStream) {
        final JSONArray propNames = new JSONArray();
        final JSONArray propBase = new JSONArray();
        final JSONArray propScale = new JSONArray();
        for(int i = 0; i < propertyNames.size(); ++i) {
            propNames.add(propertyNames.get(i));
            propBase.add(baselines.get(i));
            propScale.add(scales.get(i));
        }
        final JSONObject jsonObject = new JSONObject();
        jsonObject.put(PROPERTY_NAMES_KEY, propNames);
        jsonObject.put(PROPERTY_BASELINE_KEY, propBase);
        jsonObject.put(PROPERTY_SCALE_KEY, propScale);

        try {
            outputStream.write(jsonObject.toJSONString().getBytes());
        } catch(IOException ioException) {
            throw new GATKException("Error saving data summary json", ioException);
        }
    }

    class PropertiesTableBuilder {
        static final int initialArrayLength = 100;
        private final List<Object> properties = new ArrayList<>();
        private final List<String> propertyNames = new ArrayList<>();
        private final List<String> propertyClasses = new ArrayList<>();
        private final List<Integer> numRows = new ArrayList<>();
        private final List<Integer> maxRows = new ArrayList<>();
        private final Map<String, Float> baselines = new HashMap<>();
        private final Map<String, Float> scales = new HashMap<>();

        PropertiesTableBuilder() {}

        public PropertiesTable build() {
            final PropertiesTable propertiesTable = new PropertiesTable();
            // validate
            //   check that numRows and maxRows are the same for every property
            //   check baselines is empty or has the same keys as propertyNames (no other sensible workflow)
            //      if baselines is empty, then comptute baselines and scales
            //   resize arrays (maybe need to re-think doubler and have a generic re-sizer?)
            // add all the properties
            // return propertiesTable;
        }

        public void add(final String propertyName, final boolean value) {
            final Function<Object, Object> doubler =
                arrObj -> (Object)Arrays.copyOf((boolean[])arrObj, ((boolean[])arrObj).length * 2);
            add(propertyName, value, BOOLEAN_ARR, boolean[]::new, doubler);
        }

        public void add(final String propertyName, final boolean[] value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((boolean[][])arrObj, ((boolean[][])arrObj).length * 2);
            add(propertyName, value, BOOLEAN_MAT, boolean[][]::new, doubler);
        }

        public void add(final String propertyName, final int value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((int[])arrObj, ((int[])arrObj).length * 2);
            add(propertyName, value, INT_ARR, int[]::new, doubler);
        }

        public void add(final String propertyName, final int[] value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((int[][])arrObj, ((int[][])arrObj).length * 2);
            add(propertyName, value, INT_MAT, int[][]::new, doubler);
        }

        public void add(final String propertyName, final long value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((long[])arrObj, ((long[])arrObj).length * 2);
            add(propertyName, value, LONG_ARR, long[]::new, doubler);
        }

        public void add(final String propertyName, final long[] value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((long[][])arrObj, ((long[][])arrObj).length * 2);
            add(propertyName, value, LONG_MAT, long[][]::new, doubler);
        }

        public void add(final String propertyName, final float value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((float[])arrObj, ((float[])arrObj).length * 2);
            add(propertyName, value, FLOAT_ARR, float[]::new, doubler);
        }

        public void add(final String propertyName, final float[] value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((float[][])arrObj, ((float[][])arrObj).length * 2);
            add(propertyName, value, FLOAT_MAT, float[][]::new, doubler);
        }

        public void add(final String propertyName, final double value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((double[])arrObj, ((double[])arrObj).length * 2);
            add(propertyName, value, DOUBLE_ARR, double[]::new, doubler);
        }

        public void add(final String propertyName, final double[] value) {
            final Function<Object, Object> doubler =
                    arrObj -> (Object)Arrays.copyOf((double[][])arrObj, ((double[][])arrObj).length * 2);
            add(propertyName, value, DOUBLE_MAT, double[][]::new, doubler);
        }

        private void add(final String propertyName, final Object value, final String className,
                            final IntFunction<Object> producer, final Function<Object, Object> doubler) {
            // Maintain propertyNames in sorted order, and keep other entries in same order
            final int propertyIndex = Collections.binarySearch(propertyNames, propertyName);
            final Object propertyObj;
            if(propertyIndex < 0) {
                // first time adding property, initialize it
                propertyObj = producer.apply(initialArrayLength);
                final int insertIndex = -1 - propertyIndex;
                propertyNames.add(insertIndex, propertyName);
                properties.add(insertIndex, propertyObj);
                propertyClasses.add(insertIndex, className);
                numRows.add(insertIndex, 1);
                maxRows.add(insertIndex, initialArrayLength);
            } else {
                // append to existing property
                if(numRows.get(propertyIndex) >= maxRows.get(propertyIndex)) {
                    // array is full, reallocate with twice as much space
                    propertyObj = doubler.apply(properties.get(propertyIndex));
                    properties.set(propertyIndex, propertyObj);
                    maxRows.set(propertyIndex, maxRows.get(propertyIndex) * 2);
                } else {
                    propertyObj = properties.get(propertyIndex);
                }
                // increment number of rows
                numRows.set(propertyIndex, numRows.get(propertyIndex) + 1);
            }
            propertyObj  // need to assign value to numRows - 1
        }

        public void setBaselineAndScale(final String propertyName, final double baseline, final double scale) {
            this.baselines.put(propertyName, (float)baseline);
            this.scales.put(propertyName, (float)scale);
        }

        public void loadNormalizations(final InputStream inputStream) {
            final JSONObject jsonObject;
            try {
                jsonObject = (JSONObject) JSONValue.parseWithException(inputStream);
            } catch (IOException | ParseException ioException) {
                throw new GATKException("Unable to parse JSON from inputStream", ioException);
            }
            final JSONArray propNames = ((JSONArray) jsonObject.get(PROPERTY_NAMES_KEY));
            final JSONArray propBase = ((JSONArray) jsonObject.get(PROPERTY_BASELINE_KEY));
            final JSONArray propScale = ((JSONArray) jsonObject.get(PROPERTY_SCALE_KEY));
            for (int idx = 0; idx < propNames.size(); ++idx) {
                setBaselineAndScale(
                    (String)propNames.get(idx), getDoubleFromJSON(propBase.get(idx)), getDoubleFromJSON(propScale.get(idx))
                );
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
    }
}