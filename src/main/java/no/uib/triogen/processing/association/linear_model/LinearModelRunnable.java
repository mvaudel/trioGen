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
import org.apache.commons.math3.linear.SingularMatrixException;
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
    private final static double[] na2 = new double[]{Double.NaN, Double.NaN};
    private final static double[] na3 = new double[]{Double.NaN, Double.NaN, Double.NaN};
    private final static double[] na4 = new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN};
    private final static double[] na5 = new double[]{Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
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
        int[] hMax = new int[4];
        TreeMap<Integer, Integer>[] hHist = new TreeMap[4];

        for (int j = 0; j < 4; j++) {

            hMin[j] = Integer.MAX_VALUE;
            hMax[j] = Integer.MIN_VALUE;
            hHist[j] = new TreeMap<>();

        }

        int childMin = Integer.MAX_VALUE;
        int childMax = Integer.MIN_VALUE;
        int motherMin = Integer.MAX_VALUE;
        int motherMax = Integer.MIN_VALUE;
        int fatherMin = Integer.MAX_VALUE;
        int fatherMax = Integer.MIN_VALUE;

        TreeMap<Integer, Integer> childHist = new TreeMap<>();
        TreeMap<Integer, Integer> motherHist = new TreeMap<>();
        TreeMap<Integer, Integer> fatherHist = new TreeMap<>();

        int phenoI = 0;
        int iterationI = 0;

        for (String childId : childToParentMap.children) {

            double y = phenotypes[phenoI++];

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

                int nAltChild = h[0] + h[2];

                if (nAltChild < childMin) {

                    childMin = nAltChild;

                }
                if (nAltChild > childMax) {

                    childMax = nAltChild;

                }

                Integer frequency = childHist.get(nAltChild);

                if (frequency != null) {

                    childHist.put(nAltChild, frequency + 1);

                } else {

                    childHist.put(nAltChild, 1);

                }

                int nAltMother = h[0] + h[1];

                if (nAltMother < motherMin) {

                    motherMin = nAltMother;

                }
                if (nAltMother > motherMax) {

                    motherMax = nAltMother;

                }

                frequency = motherHist.get(nAltMother);

                if (frequency != null) {

                    motherHist.put(nAltMother, frequency + 1);

                } else {

                    motherHist.put(nAltMother, 1);

                }

                int nAltFather = h[2] + h[3];

                if (nAltFather < fatherMin) {

                    fatherMin = nAltFather;

                }
                if (nAltFather > fatherMax) {

                    fatherMax = nAltFather;

                }

                frequency = fatherHist.get(nAltFather);

                if (frequency != null) {

                    fatherHist.put(nAltFather, frequency + 1);

                } else {

                    fatherHist.put(nAltFather, 1);

                }

                hX[iterationI][0] = h[0];
                hX[iterationI][1] = h[1];
                hX[iterationI][2] = h[2];
                hX[iterationI][3] = h[3];

                cmfX[iterationI][0] = nAltChild;
                cmfX[iterationI][1] = nAltMother;
                cmfX[iterationI][2] = nAltFather;

                cmX[iterationI][0] = nAltChild;
                cmX[iterationI][1] = nAltMother;

                cfX[iterationI][0] = nAltChild;
                cfX[iterationI][1] = nAltFather;

                mfX[iterationI][0] = nAltMother;
                mfX[iterationI][1] = nAltFather;

                cX[iterationI][0] = nAltChild;

                mX[iterationI][0] = nAltMother;

                fX[iterationI][0] = nAltFather;

                phenoY[iterationI] = y;

                iterationI++;

            }
        }

        // Try to anticipate singularities
        boolean hNotSingluar = hMax[0] - hMin[0] > 0 && hMax[1] - hMin[1] > 0 && hMax[2] - hMin[2] > 0 && hMax[3] - hMin[3] > 0;
        boolean childNotSingular = childMax - childMin > 0;
        boolean motherNotSingular = motherMax - motherMin > 0;
        boolean fatherNotSingular = fatherMax - fatherMin > 0;

        if (!x0 || hNotSingluar) {

            // Build histograms
            String altHistograms = String.join("",
                    "child(",
                    getHistogramAsString(childHist),
                    ");mother(",
                    getHistogramAsString(motherHist),
                    ");father(",
                    getHistogramAsString(fatherHist),
                    ")"
            );
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
            double[] betas = na4;
            double[] betaStandardErrors = na4;
            double[] betaResiduals = na4;
            
            // Hs
            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();

            boolean regressionSuccess = false;
            if (hNotSingluar) {

                try {

                    regression.newSampleData(phenoY, hX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] hBetas;
            double[] hBetaStandardErrors;
            double hRSS;
            if (regressionSuccess) {

                hBetas = betas;
                hBetaStandardErrors = betaStandardErrors;
                double[] hBetaResiduals = betaResiduals;
                hRSS = getRSS(hBetaResiduals);

            } else {

                hBetas = na5;
                hBetaStandardErrors = na5;
                hRSS = Double.NaN;

            }

            // Child Mother Father
            regressionSuccess = false;
            if (childNotSingular && motherNotSingular && fatherNotSingular) {

                try {

                    regression.newSampleData(phenoY, cmfX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] cmfBetas;
            double[] cmfBetaStandardErrors;
            double cmfRSS;
            if (regressionSuccess) {

                cmfBetas = betas;
                cmfBetaStandardErrors = betaStandardErrors;
                double[] cmfBetaResiduals = betaResiduals;
                cmfRSS = getRSS(cmfBetaResiduals);

            } else {

                cmfBetas = na4;
                cmfBetaStandardErrors = na4;
                cmfRSS = Double.NaN;

            }

            // Child Mother
            regressionSuccess = false;
            if (childNotSingular && motherNotSingular) {

                try {

                    regression.newSampleData(phenoY, cmX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] cmBetas;
            double[] cmBetaStandardErrors;
            double cmRSS;
            if (regressionSuccess) {

                cmBetas = betas;
                cmBetaStandardErrors = betaStandardErrors;
                double[] cmBetaResiduals = betaResiduals;
                cmRSS = getRSS(cmBetaResiduals);

            } else {

                cmBetas = na3;
                cmBetaStandardErrors = na3;
                cmRSS = Double.NaN;

            }

            // Child Father
            regressionSuccess = false;
            if (childNotSingular && fatherNotSingular) {

                try {

                    regression.newSampleData(phenoY, cfX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] cfBetas;
            double[] cfBetaStandardErrors;
            double cfRSS;
            if (regressionSuccess) {

                cfBetas = betas;
                cfBetaStandardErrors = betaStandardErrors;
                double[] cfBetaResiduals = betaResiduals;
                cfRSS = getRSS(cfBetaResiduals);

            } else {

                cfBetas = na3;
                cfBetaStandardErrors = na3;
                cfRSS = Double.NaN;

            }

            // Mother Father
            regressionSuccess = false;
            if (motherNotSingular && fatherNotSingular) {

                try {

                    regression.newSampleData(phenoY, mfX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] mfBetas;
            double[] mfBetaStandardErrors;
            double mfRSS;
            if (regressionSuccess) {

                mfBetas = betas;
                mfBetaStandardErrors = betaStandardErrors;
                double[] mfBetaResiduals = betaResiduals;
                mfRSS = getRSS(mfBetaResiduals);

            } else {

                mfBetas = na3;
                mfBetaStandardErrors = na3;
                mfRSS = Double.NaN;

            }

            // Child
            regressionSuccess = false;
            if (childNotSingular) {

                try {

                    regression.newSampleData(phenoY, cX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] cBetas;
            double[] cBetaStandardErrors;
            double cRSS;
            if (regressionSuccess) {

                cBetas = betas;
                cBetaStandardErrors = betaStandardErrors;
                double[] cBetaResiduals = betaResiduals;
                cRSS = getRSS(cBetaResiduals);

            } else {

                cBetas = na2;
                cBetaStandardErrors = na2;
                cRSS = Double.NaN;

            }

            // Mother
            regressionSuccess = false;
            if (motherNotSingular) {

                try {

                    regression.newSampleData(phenoY, mX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] mBetas;
            double[] mBetaStandardErrors;
            double mRSS;
            if (regressionSuccess) {

                mBetas = betas;
                mBetaStandardErrors = betaStandardErrors;
                double[] mBetaResiduals = betaResiduals;
                mRSS = getRSS(mBetaResiduals);

            } else {

                mBetas = na2;
                mBetaStandardErrors = na2;
                mRSS = Double.NaN;

            }

            // Father
            regressionSuccess = false;
            if (fatherNotSingular) {

                try {

                    regression.newSampleData(phenoY, fX);
                    betas = regression.estimateRegressionParameters();
                    betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                    betaResiduals = regression.estimateResiduals();
                    regressionSuccess = true;

                } catch (SingularMatrixException singularMatrixException) {

                }
            }

            double[] fBetas;
            double[] fBetaStandardErrors;
            double fRSS;
            if (regressionSuccess) {

                fBetas = betas;
                fBetaStandardErrors = betaStandardErrors;
                double[] fBetaResiduals = betaResiduals;
                fRSS = getRSS(fBetaResiduals);

            } else {

                fBetas = na2;
                fBetaStandardErrors = na2;
                fRSS = Double.NaN;

            }

            // Estimate model significance
            double cmf_hP = getModelSignificance(cmfRSS, 4, hRSS, 5, nValidValues);
            double cm_cmfP = getModelSignificance(cmRSS, 3, cmfRSS, 4, nValidValues);
            double cf_cmfP = getModelSignificance(cfRSS, 3, cmfRSS, 4, nValidValues);
            double mf_cmfP = getModelSignificance(mfRSS, 3, cmfRSS, 4, nValidValues);
            double c_cmfP = getModelSignificance(cRSS, 2, cmfRSS, 4, nValidValues);
            double c_cmP = getModelSignificance(cRSS, 2, cmRSS, 3, nValidValues);
            double c_cfP = getModelSignificance(cRSS, 2, cfRSS, 3, nValidValues);
            double m_cmfP = getModelSignificance(mRSS, 2, cmfRSS, 4, nValidValues);
            double m_cmP = getModelSignificance(mRSS, 2, cmRSS, 3, nValidValues);
            double m_mfP = getModelSignificance(mRSS, 2, mfRSS, 3, nValidValues);
            double f_cmfP = getModelSignificance(fRSS, 2, cmfRSS, 4, nValidValues);
            double f_cfP = getModelSignificance(fRSS, 2, cfRSS, 3, nValidValues);
            double f_mfP = getModelSignificance(fRSS, 2, mfRSS, 3, nValidValues);

            // Estimate slope significance
            double[] hBetaP = IntStream.range(0, hBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(hBetas[j], hBetaStandardErrors[j], nValidValues - hBetas.length)
                    )
                    .toArray();
            double[] cmfBetaP = IntStream.range(0, cmfBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(cmfBetas[j], cmfBetaStandardErrors[j], nValidValues - cmfBetas.length)
                    )
                    .toArray();
            double[] cmBetaP = IntStream.range(0, cmBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(cmBetas[j], cmBetaStandardErrors[j], nValidValues - cmBetas.length)
                    )
                    .toArray();
            double[] cfBetaP = IntStream.range(0, cfBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(cfBetas[j], cfBetaStandardErrors[j], nValidValues - cfBetas.length)
                    )
                    .toArray();
            double[] mfBetaP = IntStream.range(0, mfBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(mfBetas[j], mfBetaStandardErrors[j], nValidValues - mfBetas.length)
                    )
                    .toArray();
            double[] cBetaP = IntStream.range(0, cBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(cBetas[j], cBetaStandardErrors[j], nValidValues - cBetas.length)
                    )
                    .toArray();
            double[] mBetaP = IntStream.range(0, mBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(mBetas[j], mBetaStandardErrors[j], nValidValues - mBetas.length)
                    )
                    .toArray();
            double[] fBetaP = IntStream.range(0, fBetas.length)
                    .mapToDouble(
                            j -> getBetaSignificance(fBetas[j], fBetaStandardErrors[j], nValidValues - fBetas.length)
                    )
                    .toArray();

            // Export
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder
                    .append(phenoName)
                    .append(Utils.separator)
                    .append(genotypesProvider.getVariantID())
                    .append(Utils.separator)
                    .append(nValidValues)
                    .append(Utils.separator)
                    .append(altHistograms)
                    .append(Utils.separator)
                    .append(hHistograms);
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

    /**
     * Appends the betas with se and significance to the stringBuilder with
     * preceeding separators.
     *
     * @param stringBuilder the string builder
     * @param betas estimation of the betas
     * @param betaStandardErrors standard errors of the beta estimation
     * @param betaP significance of the beta estimations
     */
    private void appendBetasAsString(
            StringBuilder stringBuilder,
            double[] betas,
            double[] betaStandardErrors,
            double[] betaP
    ) {

        for (int i = 0; i < betas.length - 1; i++) {

            stringBuilder.append(Utils.separator);
            stringBuilder.append(betas[i + 1]);
            stringBuilder.append(Utils.separator);
            stringBuilder.append(betaStandardErrors[i + 1]);
            stringBuilder.append(Utils.separator);
            stringBuilder.append(betaP[i + 1]);

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

        if (Double.isNaN(model1RSS) || Double.isNaN(model2RSS)) {

            return Double.NaN;

        }

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

        if (Double.isNaN(beta) || Double.isNaN(betaSE)) {

            return Double.NaN;

        }

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
