package no.uib.triogen.transmission.linear_model;

import java.util.ArrayList;
import java.util.stream.IntStream;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

/**
 * Runnable for the linear model association.
 *
 * @author Marc Vaudel
 */
public class LinearModelRunnable implements Runnable {

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
            double[] xValues,
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
            double[] xValues,
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

        int i = 2;
        double p = 0.0;
        double inverseAbsoluteAccuracy = 1.0;

        if (n >= 3) {

            while (p < inverseAbsoluteAccuracy) {

                long k = 1 << ++i;
                inverseAbsoluteAccuracy = FastMath.pow(10, -k);
                
                if (inverseAbsoluteAccuracy <= 0.0) {
                    
                    break;
                    
                }
                
                TDistribution tDistribution = new TDistribution(n - 2, inverseAbsoluteAccuracy);

                p = 2d * (1.0 - tDistribution.cumulativeProbability(
                        FastMath.abs(slope) / slopeSE
                ));
                
            }

        } else {

            p = Double.NaN;

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
                Long.toString(simpleRegression.getN())
        );
        outputWriter.writeLine(line);

    }
}
