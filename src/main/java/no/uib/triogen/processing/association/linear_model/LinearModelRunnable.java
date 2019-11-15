package no.uib.triogen.processing.association.linear_model;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 * Runnable for the linear model association.
 *
 * @author Marc Vaudel
 */
public class LinearModelRunnable implements Runnable {

    /**
     * Boolean indicating whether regression results should only be reported
     * when more than one value of x is available.
     */
    public static boolean x0 = false;
    /**
     * Epsilon to use for the estimation of the p-value.
     */
    private final static double[] epsilons = new double[]{1e-14, 1e-20, 1e-50, 1e-100, 1e-200};
    /**
     * Names of the hs.
     */
    private final static String[] hNames = new String[]{"h1", "h2", "h3", "h4"};
    /**
     * The variants iterator.
     */
    private final VariantIterator iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The handler for phenotypes.
     */
    private final PhenotypesHandler phenotypesHandler;
    /**
     * The output writer.
     */
    private final SimpleFileWriter outputWriter;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;

    /**
     * Constructor.
     *
     * @param iterator the variants iterator
     * @param childToParentMap the child to parent map
     * @param phenotypesHandler the phenotypes handler
     * @param outputWriter the output writer
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            PhenotypesHandler phenotypesHandler,
            SimpleFileWriter outputWriter
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.phenotypesHandler = phenotypesHandler;
        this.outputWriter = outputWriter;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider tempGenotypesProvider;
            while ((tempGenotypesProvider = iterator.next()) != null && !canceled) {

                GenotypesProvider genotypesProvider = tempGenotypesProvider;
                genotypesProvider.parse();

                ArrayList<int[]> xValuesList = new ArrayList<>(4);

                for (int i = 0; i < 4; i++) {

                    xValuesList.add(new int[childToParentMap.children.size()]);

                }

                int i = 0;

                for (String childId : childToParentMap.children) {

                    int[] hs = genotypesProvider.getH(childToParentMap, childId);

                    for (int j = 0; j < 4; j++) {

                        xValuesList.get(j)[i] = hs[j];

                    }

                    i++;

                }

                IntStream.range(0, 4)
                        .parallel()
                        .forEach(
                                hI -> runLinearModel(
                                        hI,
                                        xValuesList.get(hI),
                                        genotypesProvider
                                )
                        );
            }

        } catch (Throwable t) {

            canceled = true;
            t.printStackTrace();

        }
    }

    /**
     * Runs the linear model for the given transmission model.
     *
     * @param hI the index of the h
     * @param xValues the x values to use for regression
     * @param genotypesProvider the genotypes provider
     */
    public void runLinearModel(
            int hI,
            int[] xValues,
            GenotypesProvider genotypesProvider
    ) {

        phenotypesHandler.phenoMap.entrySet()
                .parallelStream()
                .forEach(
                        entry -> runLinearModel(
                                hI,
                                genotypesProvider,
                                xValues,
                                entry.getKey(),
                                entry.getValue()
                        )
                );
    }

    /**
     * Runs the linear model for a transmission model and a phenotype name.
     *
     * @param hI the index of the h
     * @param genotypesProvider the genotypes provider
     * @param xValues the x values to use for regression
     * @param phenoName the phenotype name
     * @param phenotypes the phenotype values
     */
    public void runLinearModel(
            int hI,
            GenotypesProvider genotypesProvider,
            int[] xValues,
            String phenoName,
            double[] phenotypes
    ) {

        // Histogram
        TreeMap<Integer, Integer> hHist = new TreeMap<>();

        // Select non-missing data points
        SimpleRegression simpleRegression = new SimpleRegression();
        int xMin = xValues[0];
        int xMax = xValues[0];

        for (int i = 0; i < xValues.length; i++) {

            int x = xValues[i];

            if (x > xMax) {

                xMax = x;

            }
            if (x < xMin) {

                xMin = x;

            }

            double y = phenotypes[i];

            if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                // Regression
                simpleRegression.addData(x, y);

                // Histogram
                Integer frequency = hHist.get(x);

                if (frequency != null) {

                    hHist.put(x, frequency + 1);

                } else {

                    hHist.put(x, 1);

                }
            }
        }

        int dx = xMax - xMin;

        // Get histogram as string
        String hHistString = hHist.entrySet().stream()
                .map(
                        entry -> String.join(
                                ":",
                                Integer.toString(entry.getKey()),
                                Integer.toString(entry.getValue())
                        )
                )
                .collect(
                        Collectors.joining(",")
                );

        // Regression results
        long n = simpleRegression.getN();
        double slope = Double.NaN;
        double slopeSE = Double.NaN;
        double p = Double.NaN;

        if (dx > 0) {

            slope = simpleRegression.getSlope();
            slopeSE = simpleRegression.getSlopeStdErr();

            long degreesOfFreedom = n - 2;

            if (degreesOfFreedom > 1 && !Double.isNaN(slope) && slopeSE > 0.0) {

                double x = slope / slopeSE;

                for (double epsilon : epsilons) {

                    p = x != 0.0
                            ? Beta.regularizedBeta(
                                    degreesOfFreedom / (degreesOfFreedom + (x * x)),
                                    0.5 * degreesOfFreedom,
                                    0.5,
                                    epsilon)
                            : 0.5;

                    if (p > epsilon * 16) {
                        break;
                    }
                }
            }
        }

        if (dx > 0 || !x0) {

            // Export
            String line = String.join(
                    Utils.separator,
                    phenoName,
                    genotypesProvider.getVariantID(),
                    hNames[hI],
                    Double.toString(slope),
                    Double.toString(slopeSE),
                    Double.toString(simpleRegression.getRSquare()),
                    Double.toString(p),
                    hHistString,
                    Long.toString(simpleRegression.getN())
            );
            outputWriter.writeLine(line);

        }
    }
}
