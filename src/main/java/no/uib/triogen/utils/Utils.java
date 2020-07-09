package no.uib.triogen.utils;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Utilities.
 *
 * @author Marc Vaudel
 */
public class Utils {

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

}
