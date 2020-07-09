package no.uib.triogen.processing.association.linear_model;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.utils.SimpleSemaphore;
import no.uib.triogen.utils.Utils;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

/**
 * Runnable for the linear model association.
 *
 * @author Marc Vaudel
 */
public class LinearModelRunnable implements Runnable {

    /**
     * If true, matrices yielding singularities will be exported.
     */
    private static final boolean debugSingularities = true;
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
     * The maf threshold. maf is computed in parents for trios where a phenotype
     * is available and values lower than threshold are not included
     * (inclusive).
     */
    private final double mafThreshold;
    /**
     * If true, dosages will be used where possible, hard calls otherwise.
     */
    private final boolean useDosages;
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
     * Mutex to keep gz file and index synchronized.
     */
    private final SimpleSemaphore gzIndexMutex;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;
    /**
     * The variants to process.
     */
    private final VariantList variantList;

    /**
     * Constructor.
     *
     * @param iterator The variants iterator.
     * @param variantList The variants to process.
     * @param mafThreshold The maf threshold. maf is computed in parents for
     * trios where a phenotype is available and values lower than threshold are
     * not included (inclusive).
     * @param useDosages If true, dosages will be used when possible, hard calls
     * otherwise.
     * @param childToParentMap The child to parent map.
     * @param models The list of the names of the models to use.
     * @param phenotypesHandler The phenotypes handler.
     * @param outputWriter The output writer.
     * @param resultsIndex The writer for the index of the results file.
     * @param gzIndexSemaphore The semaphore to keep gz file and index
     * synchronized.
     * @param logger The logger.
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            VariantList variantList,
            double mafThreshold,
            boolean useDosages,
            ChildToParentMap childToParentMap,
            Model[] models,
            PhenotypesHandler phenotypesHandler,
            IndexedGzWriter outputWriter,
            SimpleFileWriter resultsIndex,
            SimpleSemaphore gzIndexSemaphore,
            SimpleCliLogger logger
    ) {

        this.iterator = iterator;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.models = models;
        this.mafThreshold = mafThreshold;
        this.useDosages = useDosages;
        this.phenotypesHandler = phenotypesHandler;
        this.outputWriter = outputWriter;
        this.resultsIndex = resultsIndex;
        this.gzIndexMutex = gzIndexSemaphore;
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

        // Gather valid values, check maf and transmission
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
        double sumPhenos = 0.0;
        double phenoMax = Double.NaN;
        double phenoMin = Double.NaN;

        for (int i = 0; i < childToParentMap.children.length; i++) {

            double y = phenotypes[i];

            if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                if (Double.isNaN(phenoMax) || y > phenoMax) {

                    phenoMax = y;

                }
                if (Double.isNaN(phenoMin) || y < phenoMin) {

                    phenoMin = y;

                }

                sumPhenos += y;

                String childId = childToParentMap.children[i];
                short[] h = genotypesProvider.getH(childToParentMap, childId);

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

        if (!Double.isNaN(phenoMin)) {

            if (phenoMax - phenoMin > 0) {

                double maf = ((double) nAltParent) / (4 * nValidValues);

                if (maf > 0.5) {

                    maf = 1.0 - maf;

                }

                if (maf >= mafThreshold || variantList != null && variantList.contains(genotypesProvider.getVariantID())) {

                    if (!x0 || h1_0 && h1_1 && h2_0 && h2_1 && h3_0 && h3_1 && h4_0 && h4_1) {

                        // Gather the input to use for the models, gather values for the histograms
                        double[] phenoY = new double[nValidValues];

                        double phenoMean = sumPhenos / nValidValues;
                        double rss0 = 0.0;

                        ArrayList<double[][]> modelsX = new ArrayList<>(models.length);
                        ArrayList<RegressionResult> regressionResults = new ArrayList<>(models.length);

                        for (int i = 0; i < models.length; i++) {

                            Model model = models[i];

                            double[][] x = new double[nValidValues][model.betaNames.length];
                            modelsX.add(x);

                            regressionResults.add(new RegressionResult(model));

                        }

                        TreeMap<Short, Integer>[] hHist = new TreeMap[4];

                        for (int j = 0; j < 4; j++) {

                            hHist[j] = new TreeMap<>();

                        }

                        TreeMap<Short, Integer> childHist = new TreeMap<>();
                        TreeMap<Short, Integer> motherHist = new TreeMap<>();
                        TreeMap<Short, Integer> fatherHist = new TreeMap<>();

                        for (int i = 0; i < indexes.length; i++) {

                            int index = indexes[i];

                            double y = phenotypes[index];

                            double distY = y - phenoMean;
                            rss0 += distY * distY;

                            String childId = childToParentMap.children[index];
                            String motherId = childToParentMap.getMother(childId);
                            String fatherId = childToParentMap.getFather(childId);

                            short[] h = genotypesProvider.getH(
                                    childId,
                                    motherId,
                                    fatherId
                            );

                            for (int j = 0; j < 4; j++) {

                                short hJ = h[j];

                                TreeMap<Short, Integer> hHistJ = hHist[j];
                                Integer frequency = hHistJ.get(hJ);

                                if (frequency != null) {

                                    hHistJ.put(hJ, frequency + 1);

                                } else {

                                    hHistJ.put(hJ, 1);

                                }
                            }

                            short nAltChild = (short) (h[0] + h[2]);

                            Integer frequency = childHist.get(nAltChild);

                            if (frequency != null) {

                                childHist.put(nAltChild, frequency + 1);

                            } else {

                                childHist.put(nAltChild, 1);

                            }

                            short nAltMother = (short) (h[0] + h[1]);

                            frequency = motherHist.get(nAltMother);

                            if (frequency != null) {

                                motherHist.put(nAltMother, frequency + 1);

                            } else {

                                motherHist.put(nAltMother, 1);

                            }

                            short nAltFather = (short) (h[2] + h[3]);

                            frequency = fatherHist.get(nAltFather);

                            if (frequency != null) {

                                fatherHist.put(nAltFather, frequency + 1);

                            } else {

                                fatherHist.put(nAltFather, 1);

                            }

                            for (int k = 0; k < models.length; k++) {

                                Model model = models[k];
                                double[][] x = modelsX.get(k);

                                Model.fillX(
                                        x,
                                        model,
                                        i,
                                        childId,
                                        motherId,
                                        fatherId,
                                        genotypesProvider,
                                        useDosages
                                );

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

                                // Center and scale
                                Utils.centerAndScaleColumns(x);

                                // Run regression
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
                                    regressionResult.computeModelSignificance(rss0, nValidValues);
                                    regressionResult.computeBetaSignificance(nValidValues);

                                    regressionRestultsMap.put(model.name(), regressionResult);

                                } catch (SingularMatrixException singularMatrixException) {

                                    logger.logVariant(
                                            genotypesProvider.getVariantID(),
                                            String.join(" ",
                                                    model.name(),
                                                    phenoName,
                                                    "Singularity detected",
                                                    getSingularityDebugReport(
                                                            hNotSingluar,
                                                            childNotSingular,
                                                            motherNotSingular,
                                                            fatherNotSingular
                                                    )
                                            )
                                    );

                                    if (debugSingularities) {

                                        writeSingularityDebugReport(
                                                genotypesProvider.getVariantID(),
                                                phenoName,
                                                model.name(),
                                                x,
                                                phenoY
                                        );
                                    }
                                }
                            } else {

                                logger.logVariant(
                                        genotypesProvider.getVariantID(),
                                        String.join(" ",
                                                model.name(),
                                                phenoName,
                                                "Singularity anticipated",
                                                getSingularityDebugReport(
                                                        hNotSingluar,
                                                        childNotSingular,
                                                        motherNotSingular,
                                                        fatherNotSingular
                                                )
                                        )
                                );
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
                        String genotyped = genotypesProvider.genotyped() ? "1" : "0";

                        StringBuilder stringBuilder = new StringBuilder();
                        stringBuilder
                                .append(phenoName)
                                .append(IoUtils.SEPARATOR)
                                .append(genotypesProvider.getContig())
                                .append(IoUtils.SEPARATOR)
                                .append(genotypesProvider.getBp())
                                .append(IoUtils.SEPARATOR)
                                .append(genotypesProvider.getVariantID())
                                .append(IoUtils.SEPARATOR)
                                .append(genotypesProvider.getRef())
                                .append(IoUtils.SEPARATOR)
                                .append(genotypesProvider.getAlt())
                                .append(IoUtils.SEPARATOR)
                                .append(genotyped)
                                .append(IoUtils.SEPARATOR)
                                .append(nValidValues)
                                .append(IoUtils.SEPARATOR)
                                .append(altHistograms)
                                .append(IoUtils.SEPARATOR)
                                .append(hHistograms);

                        regressionResults
                                .forEach(
                                        regressionResult -> regressionResult.appendResults(stringBuilder)
                                );

                        String line = stringBuilder
                                .append(IoUtils.LINE_SEPARATOR)
                                .toString();

                        gzIndexMutex.acquire();

                        IndexedGzCoordinates coordinates = outputWriter.append(line);

                        resultsIndex.writeLine(
                                genotypesProvider.getContig(),
                                Integer.toString(genotypesProvider.getBp()),
                                genotypesProvider.getVariantID(),
                                phenoName,
                                Integer.toString(coordinates.compressedLength),
                                Integer.toString(coordinates.uncompressedLength)
                        );

                        gzIndexMutex.release();

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
            } else {

                logger.logVariant(
                        genotypesProvider.getVariantID(),
                        String.join(" ", phenoName, "has only one value.")
                );
            }
        } else {

            logger.logVariant(
                    genotypesProvider.getVariantID(),
                    String.join(" ", phenoName, "is NaN.")
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
            TreeMap<Short, Integer> hHist
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

    /**
     * Returns a summary of the anticipated singularities.
     *
     * @param hNotSingluar Boolean indicating whether h genotypes are expected
     * yield a singularity.
     * @param childNotSingular Boolean indicating whether child genotypes are
     * expected yield a singularity.
     * @param motherNotSingular Boolean indicating whether mother genotypes are
     * expected yield a singularity.
     * @param fatherNotSingular Boolean indicating whether father genotypes are
     * expected yield a singularity.
     *
     * @return Returns a summary of the anticipated singularities.
     */
    private String getSingularityDebugReport(
            boolean hNotSingluar,
            boolean childNotSingular,
            boolean motherNotSingular,
            boolean fatherNotSingular
    ) {

        StringBuilder report = new StringBuilder();
        report.append("(");
        if (!hNotSingluar) {
            report.append("h singular");
        }
        if (!childNotSingular) {
            if (report.length() > 1) {
                report.append(", ");
            }
            report.append("child singular");
        }
        if (!motherNotSingular) {
            if (report.length() > 1) {
                report.append(", ");
            }
            report.append("mother singular");
        }
        if (!fatherNotSingular) {
            if (report.length() > 1) {
                report.append(", ");
            }
            report.append("father singular");
        }
        if (report.length() == 1) {
            report.append("no singularity anticipated");
        }
        report.append(")");

        return report.toString();

    }

    private void writeSingularityDebugReport(
            String variantId,
            String phenoName,
            String modelName,
            double[][] x,
            double[] y
    ) {

        String debugFileName = String.join("_", variantId, phenoName, modelName, "debug_singularity_x");
        File singularityDebugFile = new File(outputWriter.getFile().getParentFile(), debugFileName);

        try ( SimpleFileWriter writer = new SimpleFileWriter(singularityDebugFile, true)) {

            for (int xi = 0; xi < x.length; xi++) {

                StringBuilder sb = new StringBuilder();

                for (int xj = 0; xj < x[0].length; xj++) {

                    if (xj > 0) {

                        sb.append("\t");

                    }

                    sb.append(x[xi][xj]);

                }

                writer.writeLine(sb.toString());

            }
        }

        debugFileName = String.join("_", variantId, phenoName, modelName, "debug_singularity_y");
        singularityDebugFile = new File(outputWriter.getFile().getParentFile(), debugFileName);

        try ( SimpleFileWriter writer = new SimpleFileWriter(singularityDebugFile, true)) {

            for (int yi = 0; yi < y.length; yi++) {

                writer.writeLine(Double.toString(y[yi]));

            }
        }
    }
}
