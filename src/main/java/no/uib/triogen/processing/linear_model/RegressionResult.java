package no.uib.triogen.processing.linear_model;

import java.util.Arrays;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.utils.Utils;
import static no.uib.triogen.utils.Utils.pEpsilons;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.special.Beta;

/**
 * Placeholder for the results of a regression.
 *
 * @author Marc Vaudel
 */
public class RegressionResult {

    /**
     * The model of this regression.
     */
    public final Model model;
    /**
     * Estimates for the slopes.
     */
    public double[] beta;
    /**
     * Standard error of the slope estimate.
     */
    public double[] betaStandardError;
    /**
     * Residuals.
     */
    public double[] betaResiduals;
    /**
     * Estimates for the significance of the slopes.
     */
    public double[] betaSignificance;
    /**
     * Sum of squared residuals.
     */
    public double rss = Double.NaN;
    /**
     * Model significance relatively to an intercept.
     */
    public double modelSignificance;
    /**
     * Variance explained.
     */
    public double varianceExplained;
    /**
     * Model significance relatively to included parent models.
     */
    public final double[] modelRelativeSignificance;
    /**
     * The number of values used in the model.
     */
    public final int nValues;

    /**
     * Constructor.
     *
     * @param model The model used in this regression.
     * @param nValues The number of samples used in this regression.
     */
    public RegressionResult(
            Model model,
            int nValues
    ) {

        this.model = model;
        this.nValues = nValues;

        switch (model.betaNames.length) {

            case 1:
                beta = Utils.na1;
                betaStandardError = Utils.na1;
                betaResiduals = Utils.na1;
                betaSignificance = Utils.na1;
                break;

            case 2:
                beta = Utils.na2;
                betaStandardError = Utils.na2;
                betaResiduals = Utils.na2;
                betaSignificance = Utils.na2;
                break;

            case 3:
                beta = Utils.na3;
                betaStandardError = Utils.na3;
                betaResiduals = Utils.na3;
                betaSignificance = Utils.na3;
                break;

            case 4:
                beta = Utils.na4;
                betaStandardError = Utils.na4;
                betaResiduals = Utils.na4;
                betaSignificance = Utils.na4;
                break;

            default:
                throw new IllegalArgumentException("Only 4 betas implemented.");

        }

        modelRelativeSignificance = new double[model.includedParentModels.length];
        Arrays.fill(modelRelativeSignificance, Double.NaN);

    }

    /**
     * Sets the residual sum of squares for the given residuals.
     */
    public void computeRSS() {

        rss = Arrays.stream(betaResiduals)
                .map(
                        x -> x * x
                )
                .sum();

    }

    /**
     * Computes the significance of the model relative to an intercept.
     *
     * @param rss0 The residual sum of squares for the intercept.
     */
    public void computeModelSignificance(
            double rss0
    ) {
        
        varianceExplained = (rss0 - rss) / rss0;

        modelSignificance = getModelSignificance(
                rss0,
                1,
                rss,
                model.betaNames.length + 1
        );

    }

    /**
     * Computes the significance the beta estimates.
     */
    public void computeBetaSignificance() {

        int degreesOfFreedom = nValues - beta.length - 1;

        if (degreesOfFreedom > 1) {

            betaSignificance = Arrays.copyOf(betaSignificance, betaSignificance.length);

            for (int i = 0; i < beta.length; i++) {

                double betaEstimate = beta[i];
                double betaSE = betaStandardError[i];

                if (!Double.isNaN(betaEstimate) && !Double.isNaN(betaSE) && betaSE > 0.0) {

                    double x = betaEstimate / betaSE;
                    double p = Double.NaN;

                    for (double epsilon : pEpsilons) {

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

                    betaSignificance[i] = p;

                }
            }
        }
    }

    /**
     * Computes the model significance against available parent models.
     *
     * @param regressionRestultsMap Map of the regression results.
     */
    public void computeModelSignificance(
            HashMap<String, RegressionResult> regressionRestultsMap
    ) {

        for (int i = 0; i < model.includedParentModels.length; i++) {

            String model2 = model.includedParentModels[i];

            if (regressionRestultsMap.containsKey(model2)) {

                RegressionResult regressionResult2 = regressionRestultsMap.get(model2);
                modelRelativeSignificance[i] = getModelSignificance(
                        rss,
                        model.betaNames.length + 1,
                        regressionResult2.rss,
                        regressionResult2.model.betaNames.length + 1
                );
                
            }
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
     *
     * @return the significance
     */
    private double getModelSignificance(
            double model1RSS,
            int p1,
            double model2RSS,
            int p2
    ) {

        if (Double.isNaN(model1RSS) || Double.isNaN(model2RSS)) {

            return Double.NaN;

        }

        double numeratorDegreesOfFreedom = p2 - p1;
        double denominatorDegreesOfFreedom = nValues - p2;
        double x = ((model1RSS - model2RSS) / numeratorDegreesOfFreedom) / (model2RSS / denominatorDegreesOfFreedom);

        FDistribution fDistribution = new FDistribution(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom);

        return 1.0 - fDistribution.cumulativeProbability(x);

    }

    /**
     * Appends the regression results to the given string builder. Model
     * significance followed by betas, beta standard errors, and beta
     * significance. All preceeded by a separator.
     *
     * @param stringBuilder the string builder to write to
     */
    public void appendResults(
            StringBuilder stringBuilder
    ) {

        stringBuilder
                .append(IoUtils.SEPARATOR)
                .append(varianceExplained)
                .append(IoUtils.SEPARATOR)
                .append(modelSignificance);

        Arrays.stream(modelRelativeSignificance)
                .forEach(
                        value -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(value)
                );
        Arrays.stream(beta)
                .forEach(
                        value -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(value)
                );
        Arrays.stream(betaStandardError)
                .forEach(
                        value -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(value)
                );
        Arrays.stream(betaSignificance)
                .forEach(
                        value -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(value)
                );

    }
}
