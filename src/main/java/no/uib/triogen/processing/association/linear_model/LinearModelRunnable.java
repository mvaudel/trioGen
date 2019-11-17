package no.uib.triogen.processing.association.linear_model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
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

                phenotypesHandler.phenoMap.entrySet()
                        .parallelStream()
                        .forEach(
                                entry -> runLinearModel(
                                        genotypesProvider,
                                        entry.getKey(),
                                        entry.getValue()
                                )
                        );
            }

        } catch (Throwable t) {

            canceled = true;
            t.printStackTrace();

        }
    }

    /**
     * Runs the linear models for a phenotype name.
     *
     * @param genotypesProvider the genotypes provider
     * @param phenoName the phenotype name
     * @param phenotypes the phenotype values
     */
    private void runLinearModel(
            GenotypesProvider genotypesProvider,
            String phenoName,
            double[] phenotypes
    ) {

        // Count the number of pheno values that can be used in the regression
        int nValidValues = (int) Arrays.stream(phenotypes)
                .filter(
                        y -> !Double.isNaN(y) && !Double.isInfinite(y)
                )
                .count();

        // Gather the input to use for the models
        double[] phenoY = new double[nValidValues];
        double[][] hX = new double[nValidValues][4];
        double[][] motherX = new double[nValidValues][2];
        double[][] fatherX = new double[nValidValues][2];
        double[][] childX = new double[nValidValues][2];
        int[] hMin = new int[4];
        Arrays.fill(hMin, Integer.MAX_VALUE);
        int[] hMax = new int[4];
        Arrays.fill(hMax, Integer.MIN_VALUE);
        TreeMap<Integer, Integer>[] hHist = new TreeMap[4];

        int i = 0;

        for (String childId : childToParentMap.children) {

            double y = phenotypes[i];

            if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                int[] h = genotypesProvider.getH(childToParentMap, childId);

                for (int j = 0; j < 4; j++) {

                    int hJ = h[j];

                    if (hJ < hMin[j]) {

                        hMin[j] = hJ;

                    }
                    if (hJ > hMax[j]) {

                        hMax[j] = hJ;

                    }

                    TreeMap<Integer, Integer> hHistJ = hHist[j];
                    Integer frequency = hHistJ.get(hJ);

                    if (frequency != null) {

                        hHistJ.put(hJ, frequency + 1);

                    } else {

                        hHistJ.put(hJ, 1);

                    }
                }

                hX[i][0] = h[0];
                hX[i][1] = h[2];
                hX[i][2] = h[3];
                hX[i][3] = h[4];

                motherX[i][0] = h[1] + h[2];
                motherX[i][1] = h[2] + h[3];

                fatherX[i][0] = h[1] + h[2];
                fatherX[i][1] = h[0] + h[1];

                childX[i][0] = h[0] + h[1];
                childX[i][1] = h[2] + h[3];

                i++;

            }
        }

        if (!x0 || hMax[0] - hMin[0] > 0 && hMax[1] - hMin[1] > 0 && hMax[2] - hMin[2] > 0 && hMax[3] - hMin[3] > 0) {

            // Build histograms
            String hHistograms = String.join("",
                    "h1(",
                    getHistogramAsString(hHist[0]),
                    ");h2(",
                    getHistogramAsString(hHist[1]),
                    ");h3(",
                    getHistogramAsString(hHist[2]),
                    ");h4(",
                    getHistogramAsString(hHist[3]),
                    ")"
            );

            // Run the regressions
            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();

            regression.newSampleData(phenoY, hX);
            double[] hBeta = regression.estimateRegressionParameters();
            double[] hBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] hBetaResiduals = regression.estimateResiduals();
            double hRSS = getRSS(hBetaResiduals);

            regression.newSampleData(phenoY, motherX);
            double[] motherBeta = regression.estimateRegressionParameters();
            double[] motherBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] motherBetaResiduals = regression.estimateResiduals();
            double motherRSS = getRSS(motherBetaResiduals);

            regression.newSampleData(phenoY, fatherX);
            double[] fatherBeta = regression.estimateRegressionParameters();
            double[] fatherBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] fatherBetaResiduals = regression.estimateResiduals();
            double fatherRSS = getRSS(fatherBetaResiduals);

            regression.newSampleData(phenoY, childX);
            double[] childBeta = regression.estimateRegressionParameters();
            double[] childBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] childBetaResiduals = regression.estimateResiduals();
            double childRSS = getRSS(childBetaResiduals);
            
            // Estimate model significance
            double significanceMother = getModelSignificance(motherRSS, 2, hRSS, 4, nValidValues);
            double significanceFather = getModelSignificance(fatherRSS, 2, hRSS, 4, nValidValues);
            double significanceChild = getModelSignificance(childRSS, 2, hRSS, 4, nValidValues);

            
            // Estimate slope significance
            double[] pHBeta = IntStream.range(0, 4)
                    .mapToDouble(
                            j -> getBetaSignificance(hBeta[j], hBetaStandardErrors[j], nValidValues - 4)
                    )
                    .toArray();
            double[] pMotherBeta = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(motherBeta[j], motherBetaStandardErrors[j], nValidValues - 2)
                    )
                    .toArray();
            double[] pFatherBeta = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(fatherBeta[j], fatherBetaStandardErrors[j], nValidValues - 2)
                    )
                    .toArray();
            double[] pChildBeta = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(childBeta[j], childBetaResiduals[j], nValidValues - 2)
                    )
                    .toArray();

            // Export
            String line = String.join(
                    Utils.separator,
                    phenoName,
                    genotypesProvider.getVariantID(),
                    hHistograms,
                    Integer.toString(nValidValues),
                    Double.toString(significanceMother),
                    Double.toString(significanceFather),
                    Double.toString(significanceChild),
                    Double.toString(hBeta[0]),
                    Double.toString(hBetaStandardErrors[0]),
                    Double.toString(pHBeta[0]),
                    Double.toString(hBeta[1]),
                    Double.toString(hBetaStandardErrors[1]),
                    Double.toString(pHBeta[1]),
                    Double.toString(hBeta[2]),
                    Double.toString(hBetaStandardErrors[2]),
                    Double.toString(pHBeta[2]),
                    Double.toString(hBeta[3]),
                    Double.toString(hBetaStandardErrors[3]),
                    Double.toString(pHBeta[3]),
                    Double.toString(motherBeta[0]),
                    Double.toString(motherBetaStandardErrors[0]),
                    Double.toString(pMotherBeta[0]),
                    Double.toString(motherBeta[1]),
                    Double.toString(motherBetaStandardErrors[1]),
                    Double.toString(pMotherBeta[1]),
                    Double.toString(fatherBeta[0]),
                    Double.toString(fatherBetaStandardErrors[0]),
                    Double.toString(pFatherBeta[0]),
                    Double.toString(fatherBeta[1]),
                    Double.toString(fatherBetaStandardErrors[1]),
                    Double.toString(pFatherBeta[1]),
                    Double.toString(childBeta[0]),
                    Double.toString(childBetaStandardErrors[0]),
                    Double.toString(pChildBeta[0]),
                    Double.toString(childBeta[1]),
                    Double.toString(childBetaStandardErrors[1]),
                    Double.toString(pChildBeta[1])
            );
            outputWriter.writeLine(line);

        }
    }

    /**
     * Returns the significance of increasing the complexity of a simple model,
     * model1, with p1 parameters to a more complex model, model2, with p2
     * parameters.
     *
     * @param model1RSS the residual sum of squares for model 1
     * @param p1 the number of parameters of model 1
     * @param model2RSS the residual sum of squares for model 2
     * @param p2 the number of parameters of model 2
     * @param n the number of values
     *
     * @return the significance
     */
    private double getModelSignificance(
            double model1RSS,
            int p1,
            double model2RSS,
            int p2,
            int n
    ) {

        double numeratorDegreesOfFreedom = p2 - p1;
        double denominatorDegreesOfFreedom = n - p2;
        double x = ((model1RSS - model2RSS) / numeratorDegreesOfFreedom) / (model2RSS / denominatorDegreesOfFreedom);

        FDistribution fDistribution = new FDistribution(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom);

        return fDistribution.cumulativeProbability(x);

    }

    /**
     * Returns the residual sum of squares for the given residuals.
     *
     * @param residuals the residuals
     *
     * @return the residual sum of squares
     */
    private double getRSS(
            double[] residuals
    ) {

        return Arrays.stream(residuals)
                .map(
                        x -> x * x
                )
                .sum();

    }

    /**
     * Returns the significance of a beta.
     *
     * @param beta the estimate of beta
     * @param betaSE the standard error of the estimate
     * @param degreesOfFreedom the number of degrees of freedom
     *
     * @return the significance of a beta
     */
    private double getBetaSignificance(
            double beta,
            double betaSE,
            int degreesOfFreedom
    ) {

        double p = Double.NaN;

        if (degreesOfFreedom > 1 && !Double.isNaN(beta) && betaSE > 0.0) {

            double x = beta / betaSE;

            for (double epsilon : epsilons) {

                p = x != 0.0
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

        return p;

    }

    /**
     * Returns the histogram of h as a string.
     *
     * @param hHist the histogram of h
     *
     * @return a string for the histogram
     */
    private String getHistogramAsString(
            TreeMap<Integer, Integer> hHist
    ) {

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
}
