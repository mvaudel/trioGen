package no.uib.triogen.processing.association.linear_model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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

                ArrayList<double[]> xValuesList = new ArrayList<>(4);

                for (int i = 0; i < 4; i++) {

                    xValuesList.add(new double[childToParentMap.children.size()]);

                }

                int i = 0;

                for (String childId : childToParentMap.children) {

                    double[] hs = genotypesProvider.getH(childToParentMap, childId);

                    for (int j = 0; j < 4; j++) {

                        double xValue = hs[j];

                        xValuesList.get(j)[i] = xValue;

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
            double[] xValues,
            GenotypesProvider genotypesProvider
    ) {

        String hHistString = getHHist(xValues);

        phenotypesHandler.phenoMap.entrySet()
                .parallelStream()
                .forEach(
                        entry -> runLinearModel(
                                hI,
                                genotypesProvider,
                                xValues,
                                hHistString,
                                entry.getKey(),
                                entry.getValue()
                        )
                );
    }
    
    /**
     * Returns the histogram of h as string. h is here estimated by rounding x. The result of 123 h=0 and 32 h=1 is in the form 0:123,1:32.
     * 
     * @param xValues the values of h
     * 
     * @return 
     */
    public String getHHist(double[] xValues) {
        
        TreeMap<Integer, Integer> hHist = new TreeMap<>();

        for (double xValue : xValues) {

            int h = (int) Math.round(xValue);

            Integer frequency = hHist.get(h);

            if (frequency != null) {

                hHist.put(h, frequency + 1);

            } else {

                hHist.put(h, 1);

            }
        }

        return hHist.entrySet().stream()
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
    }

    /**
     * Runs the linear model for a transmission model and a phenotype name.
     *
     * @param hI the index of the h
     * @param genotypesProvider the genotypes provider
     * @param xValues the x values to use for regression
     * @param hHistString the frequencies of the h as string
     * @param phenoName the phenotype name
     * @param phenotypes the phenotype values
     */
    public void runLinearModel(
            int hI,
            GenotypesProvider genotypesProvider,
            double[] xValues,
            String hHistString,
            String phenoName,
            double[] phenotypes
    ) {

        // Run regression
        SimpleRegression simpleRegression = new SimpleRegression();

        for (int i = 0; i < xValues.length; i++) {

            double x = xValues[i];

            if (!Double.isNaN(x) && !Double.isInfinite(x)) {

                double y = phenotypes[i];

                if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                    simpleRegression.addData(x, y);

                }
            }
        }

        // Estimate significance
        long n = simpleRegression.getN();
        double slope = simpleRegression.getSlope();
        double slopeSE = simpleRegression.getSlopeStdErr();

        long degreesOfFreedom = n - 2;

        double p = Double.NaN;

        if (degreesOfFreedom > 1 && slopeSE > 0.0) {

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
