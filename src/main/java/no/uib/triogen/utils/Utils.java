package no.uib.triogen.utils;

import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Beta;

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
     * Epsilon to use for the estimation of the p-value.
     */
    public final static double[] pEpsilons = new double[]{1e-14, 1e-20, 1e-50, 1e-100, 1e-200};
    /**
     * Cache for a standard normal distribution.
     */
    private final static NormalDistribution normalDistribution = new NormalDistribution(0, 1);

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

    /**
     * Computes the significance the beta estimates.
     *
     * @param nSamples The number of samples that went into the analysis.
     * @param nVariables The number of variables in the model excluding
     * intercept.
     * @param beta The effect size estimate.
     * @param se The standard error estimate.
     *
     * @return The significance the beta estimates.
     */
    public double computeBetaSignificance(
            int nSamples,
            int nVariables,
            double beta,
            double se
    ) {

        if (!Double.isNaN(beta) && !Double.isNaN(se) && se > 0.0) {

            int degreesOfFreedom = nSamples - nVariables - 1;

            if (degreesOfFreedom > 1) {

                double x = beta / se;

                for (double epsilon : pEpsilons) {

                    double p = x != 0.0
                            ? Beta.regularizedBeta(
                                    degreesOfFreedom / (degreesOfFreedom + (x * x)),
                                    0.5 * degreesOfFreedom,
                                    0.5,
                                    epsilon)
                            : 0.5;

                    if (p > epsilon * 16) {

                        return p;

                    }
                }
            }
        }

        return Double.NaN;

    }

    /**
     * Computes the significance the beta estimates.
     *
     * @param beta The effect size estimate.
     * @param se The standard error estimate.
     *
     * @return The significance the beta estimates.
     */
    public static double computeBetaSignificance(
            double beta,
            double se
    ) {
        
        return 2* normalDistribution.cumulativeProbability(-Math.abs(beta/se));
        
    }

}
