package no.uib.triogen.processing.association.linear_model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.gz.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.gz.IndexedGzWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.Logger;
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
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included.
     */
    private final double mafThreshold;
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
    private final IndexedGzWriter outputWriter;
    /**
     * Writer for the index of the results file.
     */
    private final SimpleFileWriter resultsIndex;
    /**
     * The logger.
     */
    private final Logger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;

    /**
     * Constructor.
     *
     * @param iterator the variants iterator
     * @param mafThreshold the maf threshold
     * @param childToParentMap the child to parent map
     * @param models the list of the names of the models to use
     * @param phenotypesHandler the phenotypes handler
     * @param outputWriter the output writer
     * @param resultsIndex writer for the index of the results file
     * @param logger the logger
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            double mafThreshold,
            ChildToParentMap childToParentMap,
            Model[] models,
            PhenotypesHandler phenotypesHandler,
            IndexedGzWriter outputWriter,
            SimpleFileWriter resultsIndex,
            Logger logger
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.models = models;
        this.mafThreshold = mafThreshold;
        this.phenotypesHandler = phenotypesHandler;
        this.outputWriter = outputWriter;
        this.resultsIndex = resultsIndex;
        this.logger = logger;

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
                                        entry.getValue(),
                                        phenotypesHandler.nValidValuesMap.get(entry.getKey())
                                )
                        );
            }

        } catch (Throwable t) {

            canceled = true;

            logger.logError(
                    Arrays.stream(t.getStackTrace())
                            .map(
                                    element -> element.toString()
                            )
                            .collect(Collectors.joining(" "))
            );

            t.printStackTrace();

        }
    }

    /**
     * Runs the linear models for a phenotype name.
     *
     * @param genotypesProvider the genotypes provider
     * @param phenoName the phenotype name
     * @param phenotypes the phenotype values
     * @param nValidValues the number of phenotypes that are finite and not NA
     */
    private void runLinearModel(
            GenotypesProvider genotypesProvider,
            String phenoName,
            double[] phenotypes,
            int nValidValues
    ) {

        // Gether valid values, check maf and transmission
        int[] indexes = new int[nValidValues];

        int cpt = 0;
        int nAltParent = 0;
        boolean h1_0 = false;
        boolean h1_1 = false;
        boolean h2_0 = false;
        boolean h2_1 = false;
        boolean h3_0 = false;
        boolean h3_1 = false;
        boolean h4_0 = false;
        boolean h4_1 = false;

        for (int i = 0; i < childToParentMap.children.length; i++) {

            double y = phenotypes[i];

            if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                String childId = childToParentMap.children[i];
                int[] h = genotypesProvider.getH(childToParentMap, childId);

                h1_0 = h1_0 || h[0] == 0;
                h1_1 = h1_1 || h[0] == 1;
                h2_0 = h2_0 || h[1] == 0;
                h2_1 = h2_1 || h[1] == 1;
                h3_0 = h3_0 || h[2] == 0;
                h3_1 = h3_1 || h[2] == 1;
                h4_0 = h4_0 || h[3] == 0;
                h4_1 = h4_1 || h[3] == 1;

                nAltParent += h[0] + h[1] + h[2] + h[3];

                indexes[cpt++] = i;

            }
        }

        double maf = ((double) nAltParent) / (4 * nValidValues);

        if (maf > mafThreshold) {

            if (!x0 || h1_0 && h1_1 && h2_0 && h2_1 && h3_0 && h3_1 && h4_0 && h4_1) {

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

                TreeMap<Integer, Integer>[] hHist = new TreeMap[4];

                for (int j = 0; j < 4; j++) {

                    hHist[j] = new TreeMap<>();

                }

                TreeMap<Integer, Integer> childHist = new TreeMap<>();
                TreeMap<Integer, Integer> motherHist = new TreeMap<>();
                TreeMap<Integer, Integer> fatherHist = new TreeMap<>();

                for (int i = 0; i < indexes.length; i++) {

                    int index = indexes[i];

                    double y = phenotypes[index];
                    String childId = childToParentMap.children[index];

                    int[] h = genotypesProvider.getH(childToParentMap, childId);

                    for (int j = 0; j < 4; j++) {

                        int hJ = h[j];

                        TreeMap<Integer, Integer> hHistJ = hHist[j];
                        Integer frequency = hHistJ.get(hJ);

                        if (frequency != null) {

                            hHistJ.put(hJ, frequency + 1);

                        } else {

                            hHistJ.put(hJ, 1);

                        }
                    }

                    int nAltChild = h[0] + h[2];

                    Integer frequency = childHist.get(nAltChild);

                    if (frequency != null) {

                        childHist.put(nAltChild, frequency + 1);

                    } else {

                        childHist.put(nAltChild, 1);

                    }

                    int nAltMother = h[0] + h[1];

                    frequency = motherHist.get(nAltMother);

                    if (frequency != null) {

                        motherHist.put(nAltMother, frequency + 1);

                    } else {

                        motherHist.put(nAltMother, 1);

                    }

                    int nAltFather = h[2] + h[3];

                    frequency = fatherHist.get(nAltFather);

                    if (frequency != null) {

                        fatherHist.put(nAltFather, frequency + 1);

                    } else {

                        fatherHist.put(nAltFather, 1);

                    }

                    for (int k = 0; k < models.length; k++) {

                        Model model = models[k];
                        double[][] x = modelsX.get(k);

                        Model.fillX(x, model, i, h, nAltChild, nAltMother, nAltFather);

                    }

                    phenoY[i] = y;

                }

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

                // Anticipate singularities
                boolean hNotSingluar = h1_0 && h1_1 && h2_0 && h2_1 && h3_0 && h3_1 && h4_0 && h4_1;
                boolean childNotSingular = h2_0 && h2_1 && h3_0 && h3_1;
                boolean motherNotSingular = h1_0 && h1_1 && h2_0 && h2_1;
                boolean fatherNotSingular = h3_0 && h3_1 && h4_0 && h4_1;

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

                String line = stringBuilder
                        .append(IoUtils.lineSeparator)
                        .toString();

                IndexedGzCoordinates coordinates = outputWriter.append(line);

                resultsIndex.writeLine(
                        genotypesProvider.getVariantID(),
                        phenoName,
                        Integer.toString(coordinates.compressedLength),
                        Integer.toString(coordinates.uncompressedLength)
                );

            } else {

                logger.logVariant(
                        genotypesProvider.getVariantID(),
                        "no transmission"
                );
            }
        } else {

            logger.logVariant(
                    genotypesProvider.getVariantID(),
                    "low maf"
            );
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
