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
        double[][] cmfX = new double[nValidValues][3];
        double[][] cfX = new double[nValidValues][2];
        double[][] cmX = new double[nValidValues][2];
        double[][] mfX = new double[nValidValues][2];
        double[][] cX = new double[nValidValues][1];
        double[][] mX = new double[nValidValues][1];
        double[][] fX = new double[nValidValues][1];
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

                cmfX[i][0] = h[1] + h[2];
                cmfX[i][1] = h[0] + h[1];
                cmfX[i][2] = h[2] + h[3];

                cfX[i][0] = h[1] + h[2];
                cfX[i][1] = h[2] + h[3];

                cmX[i][0] = h[1] + h[2];
                cmX[i][1] = h[0] + h[1];

                mfX[i][0] = h[0] + h[1];
                mfX[i][1] = h[2] + h[3];

                cX[i][0] = h[1] + h[2];

                mX[i][0] = h[0] + h[1];

                fX[i][0] = h[2] + h[3];

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
            double[] hBetas = regression.estimateRegressionParameters();
            double[] hBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] hBetaResiduals = regression.estimateResiduals();
            double hRSS = getRSS(hBetaResiduals);

            regression.newSampleData(phenoY, cmfX);
            double[] cmfBetas = regression.estimateRegressionParameters();
            double[] cmfBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] cmfBetaResiduals = regression.estimateResiduals();
            double cmfRSS = getRSS(cmfBetaResiduals);

            regression.newSampleData(phenoY, cmX);
            double[] cmBetas = regression.estimateRegressionParameters();
            double[] cmBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] cmBetaResiduals = regression.estimateResiduals();
            double cmRSS = getRSS(cmBetaResiduals);

            regression.newSampleData(phenoY, cfX);
            double[] cfBetas = regression.estimateRegressionParameters();
            double[] cfBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] cfBetaResiduals = regression.estimateResiduals();
            double cfRSS = getRSS(cfBetaResiduals);

            regression.newSampleData(phenoY, mfX);
            double[] mfBetas = regression.estimateRegressionParameters();
            double[] mfBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] mfBetaResiduals = regression.estimateResiduals();
            double mfRSS = getRSS(mfBetaResiduals);

            regression.newSampleData(phenoY, cX);
            double[] cBetas = regression.estimateRegressionParameters();
            double[] cBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] cBetaResiduals = regression.estimateResiduals();
            double cRSS = getRSS(cBetaResiduals);

            regression.newSampleData(phenoY, mX);
            double[] mBetas = regression.estimateRegressionParameters();
            double[] mBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] mBetaResiduals = regression.estimateResiduals();
            double mRSS = getRSS(mBetaResiduals);

            regression.newSampleData(phenoY, fX);
            double[] fBetas = regression.estimateRegressionParameters();
            double[] fBetaStandardErrors = regression.estimateRegressionParametersStandardErrors();
            double[] fBetaResiduals = regression.estimateResiduals();
            double fRSS = getRSS(fBetaResiduals);

            // Estimate model significance
            double cmf_hP = getModelSignificance(cmfRSS, 3, hRSS, 4, nValidValues);
            double cm_cmfP = getModelSignificance(cmRSS, 2, cmfRSS, 3, nValidValues);
            double cf_cmfP = getModelSignificance(cfRSS, 2, cmfRSS, 3, nValidValues);
            double mf_cmfP = getModelSignificance(mfRSS, 2, cmfRSS, 3, nValidValues);
            double c_cmfP = getModelSignificance(cRSS, 1, cmfRSS, 3, nValidValues);
            double c_cmP = getModelSignificance(cRSS, 1, cmRSS, 2, nValidValues);
            double c_cfP = getModelSignificance(cRSS, 1, cfRSS, 2, nValidValues);
            double m_cmfP = getModelSignificance(mRSS, 1, cmfRSS, 3, nValidValues);
            double m_cmP = getModelSignificance(mRSS, 1, cmRSS, 2, nValidValues);
            double m_mfP = getModelSignificance(mRSS, 1, mfRSS, 2, nValidValues);
            double f_cmfP = getModelSignificance(fRSS, 1, cmfRSS, 3, nValidValues);
            double f_cfP = getModelSignificance(fRSS, 1, cfRSS, 2, nValidValues);
            double f_mfP = getModelSignificance(fRSS, 1, mfRSS, 2, nValidValues);

            // Estimate slope significance
            double[] hBetaP = IntStream.range(0, 4)
                    .mapToDouble(
                            j -> getBetaSignificance(hBetas[j], hBetaStandardErrors[j], nValidValues - 4)
                    )
                    .toArray();
            double[] cmfBetaP = IntStream.range(0, 3)
                    .mapToDouble(
                            j -> getBetaSignificance(cmfBetas[j], cmfBetaStandardErrors[j], nValidValues - 3)
                    )
                    .toArray();
            double[] cmBetaP = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(cmBetas[j], cmBetaStandardErrors[j], nValidValues - 2)
                    )
                    .toArray();
            double[] cfBetaP = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(cfBetas[j], cfBetaStandardErrors[j], nValidValues - 2)
                    )
                    .toArray();
            double[] mfBetaP = IntStream.range(0, 2)
                    .mapToDouble(
                            j -> getBetaSignificance(mfBetas[j], mfBetaStandardErrors[j], nValidValues - 2)
                    )
                    .toArray();
            double[] cBetaP = IntStream.range(0, 1)
                    .mapToDouble(
                            j -> getBetaSignificance(cBetas[j], cBetaStandardErrors[j], nValidValues - 1)
                    )
                    .toArray();
            double[] mBetaP = IntStream.range(0, 1)
                    .mapToDouble(
                            j -> getBetaSignificance(mBetas[j], cBetaStandardErrors[j], nValidValues - 1)
                    )
                    .toArray();
            double[] fBetaP = IntStream.range(0, 1)
                    .mapToDouble(
                            j -> getBetaSignificance(fBetas[j], cBetaStandardErrors[j], nValidValues - 1)
                    )
                    .toArray();

            // Export
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder
                    .append(phenoName)
                    .append(Utils.separator)
                    .append(genotypesProvider.getVariantID())
                    .append(Utils.separator)
                    .append(hHistograms)
                    .append(Utils.separator)
                    .append(nValidValues);
            appendBetasAsString(
                    stringBuilder, 
                    hBetas, 
                    hBetaStandardErrors, 
                    hBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(cmf_hP);
            appendBetasAsString(
                    stringBuilder, 
                    cmfBetas, 
                    cmfBetaStandardErrors, 
                    cmfBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(cm_cmfP);
            appendBetasAsString(
                    stringBuilder, 
                    cmBetas, 
                    cmBetaStandardErrors, 
                    cmBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(cf_cmfP);
            appendBetasAsString(
                    stringBuilder, 
                    cfBetas, 
                    cfBetaStandardErrors, 
                    cfBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(mf_cmfP);
            appendBetasAsString(
                    stringBuilder, 
                    mfBetas, 
                    mfBetaStandardErrors, 
                    mfBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(c_cmfP)
                    .append(Utils.separator)
                    .append(c_cmP)
                    .append(Utils.separator)
                    .append(c_cfP);
            appendBetasAsString(
                    stringBuilder, 
                    cBetas, 
                    cBetaStandardErrors, 
                    cBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(m_cmfP)
                    .append(Utils.separator)
                    .append(m_cmP)
                    .append(Utils.separator)
                    .append(m_mfP);
            appendBetasAsString(
                    stringBuilder, 
                    mBetas, 
                    mBetaStandardErrors, 
                    mBetaP
            );
            stringBuilder
                    .append(Utils.separator)
                    .append(f_cmfP)
                    .append(Utils.separator)
                    .append(f_cfP)
                    .append(Utils.separator)
                    .append(f_mfP);
            appendBetasAsString(
                    stringBuilder, 
                    fBetas, 
                    fBetaStandardErrors, 
                    fBetaP
            );
            String line = stringBuilder.toString();
            outputWriter.writeLine(line);

        }
    }

    private void appendBetasAsString(
            StringBuilder stringBuilder,
            double[] betas,
            double[] betaStandardErrors,
            double[] betaP
    ) {

        for (int i = 0; i < betas.length; i++) {

            stringBuilder.append(Utils.separator);
            stringBuilder.append(betas[i]);
            stringBuilder.append(Utils.separator);
            stringBuilder.append(betaStandardErrors[i]);
            stringBuilder.append(Utils.separator);
            stringBuilder.append(betaP[i]);

        }
    }

    /**
     * Returns the significance of increasing the complexity of a simple model,
     * model1, with p1 parameters to a more complex model, model2, with p2
     * parameters. p1 < p2.
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
