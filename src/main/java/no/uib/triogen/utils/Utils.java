package no.uib.triogen.utils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Utilities.
 *
 * @author Marc Vaudel
 */
public class Utils {

    /**
     * Wild card to flag contig names in file paths.
     */
    public final static String CONTIG_WILDCARD = "{contig}";
    /**
     * Placeholder for an array of NAs.
     */
    public final static double[] na1 = new double[]{Double.NaN};
    /**
     * Placeholder for an array of NAs.
     */
    public final static double[] na2 = new double[]{Double.NaN, Double.NaN};
    /**
     * Placeholder for an array of NAs.
     */
    public final static double[] na3 = new double[]{Double.NaN, Double.NaN, Double.NaN};
    /**
     * Placeholder for an array of NAs.
     */
    public final static double[] na4 = new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};

    /**
     * Simple method to merge two byte arrays.
     *
     * @param array1 First byte array.
     * @param array2 Second byte array.
     * @param len2 The length of the second array to copy
     *
     * @return A concatenation of the first and the second arrays.
     */
    public static byte[] mergeArrays(
            byte[] array1,
            byte[] array2,
            int len2
    ) {

        byte[] result = new byte[array1.length + len2];

        System.arraycopy(array1, 0, result, 0, array1.length);
        System.arraycopy(array2, 0, result, array1.length, len2);

        return result;

    }

    /**
     * Centers and scales the values in x using the mean and standard deviation.
     *
     * @param x The array to scale.
     */
    public static void centerAndScale(
            double[] x
    ) {

        int n = x.length;

        double mean = Arrays.stream(x)
                .filter(
                        xi -> !Double.isNaN(xi)
                )
                .sum() / n;

        double sd = Math.sqrt(
                Arrays.stream(x)
                        .filter(
                                xi -> !Double.isNaN(xi)
                        )
                        .map(
                                xi -> (xi - mean) * (xi - mean)
                        )
                        .sum() / n
        );

        for (int i = 0; i < n; i++) {

            if (!Double.isNaN(x[i])) {

                x[i] = (x[i] - mean) / sd;

            }
        }
    }

    /**
     * Method to estimate the median of a sorted list.
     *
     * @param input ArrayList of double
     * @return median of the input
     */
    public static double medianSorted(ArrayList<Double> input) {
        return percentileSorted(input, 0.5);
    }

    /**
     * Returns the desired percentile in a list of sorted double values. If the
     * percentile is between two values a linear interpolation is done. The list
     * must be sorted prior to submission.
     *
     * @param input the input list
     * @param percentile the desired percentile. 0.01 returns the first
     * percentile. 0.5 returns the median.
     *
     * @return the desired percentile
     */
    public static double percentileSorted(ArrayList<Double> input, double percentile) {

        if (percentile < 0 || percentile > 1) {
            throw new IllegalArgumentException(
                    "Incorrect input for percentile: "
                    + percentile + ". Input must be between 0 and 1.");
        }

        if (input == null) {
            throw new IllegalArgumentException(
                    "Attempting to estimate the percentile of a null object.");
        }

        int length = input.size();

        if (length == 0) {
            throw new IllegalArgumentException(
                    "Attempting to estimate the percentile of an empty list.");
        }

        if (length == 1) {
            return input.get(0);
        }

        double indexDouble = percentile * (length - 1);
        int index = (int) (indexDouble);
        double valueAtIndex = input.get(index);
        double rest = indexDouble - index;

        if (index == input.size() - 1 || rest == 0) {
            return valueAtIndex;
        }

        return valueAtIndex + rest * (input.get(index + 1) - valueAtIndex);
    }

}
