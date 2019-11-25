package no.uib.triogen.processing.association.linear_model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.Model;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

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
     * The variants iterator.
     */
    private final VariantIterator iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * List of models to use.
     */
    public Model[] models;
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
     * @param models the list of the names of the models to use
     * @param phenotypesHandler the phenotypes handler
     * @param outputWriter the output writer
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            Model[] models,
            PhenotypesHandler phenotypesHandler,
            SimpleFileWriter outputWriter
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.models = models;
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
        ArrayList<double[][]> modelsX = new ArrayList<>(models.length);
        ArrayList<RegressionResult> regressionResults = new ArrayList<>(models.length);

        for (int i = 0; i < models.length; i++) {

            Model model = models[i];

            double[][] x = new double[nValidValues][model.betaNames.length];
            modelsX.add(x);

            regressionResults.add(new RegressionResult(model));

        }

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

                for (int i = 0; i < models.length; i++) {

                    Model model = models[i];
                    double[][] x = modelsX.get(i);

                    Model.fillX(x, model, iterationI, h, nAltChild, nAltMother, nAltFather);

                }

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
            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
            HashMap<String, RegressionResult> regressionRestultsMap = new HashMap<>(regressionResults.size());

            for (int i = 0; i < models.length; i++) {

                Model model = models[i];

                if (Model.likelyNotSingular(
                        model,
                        hNotSingluar,
                        childNotSingular,
                        motherNotSingular,
                        fatherNotSingular
                )) {

                    double[][] x = modelsX.get(i);
                    RegressionResult regressionResult = regressionResults.get(i);

                    try {

                        regression.newSampleData(phenoY, x);
                        double[] betas = regression.estimateRegressionParameters();
                        double[] betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                        double[] betaResiduals = regression.estimateResiduals();

                        regressionResult.beta = Arrays.copyOfRange(betas, 1, betas.length);
                        regressionResult.betaStandardError = Arrays.copyOfRange(betaStandardErrors, 1, betaStandardErrors.length);
                        regressionResult.betaResiduals = Arrays.copyOfRange(betaResiduals, 1, betaResiduals.length);

                        regressionResult.computeRSS();
                        regressionResult.computeBetaSignificance(nValidValues);

                        regressionRestultsMap.put(model.name(), regressionResult);

                    } catch (SingularMatrixException singularMatrixException) {

                    }
                }
            }

            // Estimate model significance
            regressionRestultsMap.values()
                    .forEach(
                            regressionResult -> regressionResult.computeModelSignificance(
                                    regressionRestultsMap,
                                    nValidValues
                            )
                    );

            // Export
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder
                    .append(phenoName)
                    .append(IoUtils.separator)
                    .append(genotypesProvider.getVariantID())
                    .append(IoUtils.separator)
                    .append(nValidValues)
                    .append(IoUtils.separator)
                    .append(altHistograms)
                    .append(IoUtils.separator)
                    .append(hHistograms);
            
            regressionResults
                    .forEach(
                            regressionResult -> regressionResult.appendResults(stringBuilder)
                    );
            String line = stringBuilder.toString();
            outputWriter.writeLine(line);

        }
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
