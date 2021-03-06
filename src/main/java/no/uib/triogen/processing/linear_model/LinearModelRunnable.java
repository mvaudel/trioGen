package no.uib.triogen.processing.linear_model;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.io.IoUtils;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.covariates.CovariatesHandler;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.mendelian_error.MendelianErrorEstimator;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.utils.SimpleSemaphore;
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
     * The index of the bgen file to process.
     */
    private final BgenIndex bgenIndex;
    /**
     * The reader for the bgen file to process.
     */
    private final BgenFileReader bgenFileReader;
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
    private final double alleleFrequencyThreshold;
    /**
     * List of models to use.
     */
    public Model[] models;
    /**
     * The handler for phenotypes.
     */
    private final PhenotypesHandler phenotypesHandler;
    /**
     * The handler for covariates.
     */
    private final CovariatesHandler covariatesHandler;
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
     * The decompressor to use.
     */
    private final ZstdDecompressor decompressor = new ZstdDecompressor();

    /**
     * Constructor.
     *
     * @param iterator The variants iterator.
     * @param bgenIndex The index of the bgen file.
     * @param bgenFileReader The reader for the bgen file.
     * @param variantList The variants to process.
     * @param frequencyThreshold The maf threshold. maf is computed in parents
     * for trios where a phenotype is available and values lower than threshold
     * are not included (inclusive).
     * @param childToParentMap The child to parent map.
     * @param models The list of the names of the models to use.
     * @param phenotypesHandler The phenotypes handler.
     * @param covariatesHandler The covariates handler.
     * @param outputWriter The output writer.
     * @param resultsIndex The writer for the index of the results file.
     * @param gzIndexSemaphore The semaphore to keep gz file and index
     * synchronized.
     * @param logger The logger.
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            BgenIndex bgenIndex,
            BgenFileReader bgenFileReader,
            VariantList variantList,
            double frequencyThreshold,
            ChildToParentMap childToParentMap,
            Model[] models,
            PhenotypesHandler phenotypesHandler,
            CovariatesHandler covariatesHandler,
            IndexedGzWriter outputWriter,
            SimpleFileWriter resultsIndex,
            SimpleSemaphore gzIndexSemaphore,
            SimpleCliLogger logger
    ) {

        this.iterator = iterator;
        this.bgenIndex = bgenIndex;
        this.bgenFileReader = bgenFileReader;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.models = models;
        this.alleleFrequencyThreshold = frequencyThreshold;
        this.phenotypesHandler = phenotypesHandler;
        this.covariatesHandler = covariatesHandler;
        this.outputWriter = outputWriter;
        this.resultsIndex = resultsIndex;
        this.gzIndexMutex = gzIndexSemaphore;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            Integer tempIndex;
            while ((tempIndex = iterator.next()) != null && !canceled) {

                int variantIndex = tempIndex;
                VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

                if (variantInformation.alleles.length > 1) {

                    if (variantList == null || variantList.include(variantInformation.contig, variantInformation.position)) {

                        BgenVariantData variantData = bgenFileReader.getVariantData(variantIndex);
                        variantData.parse(
                                childToParentMap,
                                decompressor
                        );

                        // Get the alleles passing the frequency threshold, test all alleles if the variant is targeted
                        int[] testedAlleleIndexes = variantList == null || !variantList.contains(variantInformation.id) && !variantList.contains(variantInformation.rsid)
                                ? IntStream.range(1, variantData.getOrderedAlleles().length)
                                        .filter(
                                                alleleIndex -> variantData.getAlleleFrequency(alleleIndex) > alleleFrequencyThreshold
                                                && variantData.getAlleleFrequency(alleleIndex) < 1.0 - alleleFrequencyThreshold
                                        )
                                        .toArray() : IntStream.range(1, variantData.getOrderedAlleles().length).toArray();

                        if (testedAlleleIndexes.length > 0) {

                            int alleleI = testedAlleleIndexes[0];

                            // Estimate the prevalence of Mendelian errors, swap child alleles if >50%
                            double mendelianErrors = MendelianErrorEstimator.estimateMendelianErrorPrevalence(variantData, childToParentMap, alleleI);

                            if (!Double.isNaN(mendelianErrors) && mendelianErrors > 0.5) {

                                variantData.swapChildrenAlleles();

                            }

                            // Run linear model
                            phenotypesHandler.phenoMap.entrySet()
                                    .parallelStream()
                                    .forEach(
                                            entry -> runLinearModel(
                                                    variantIndex,
                                                    testedAlleleIndexes,
                                                    variantData,
                                                    entry.getKey(),
                                                    entry.getValue()
                                            )
                                    );
                        }
                    }
                }
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
     * @param variantIndex The index of this variant in the bgen file.
     * @param testedAlleleIndexes The indexes of the alleles to test.
     * @param variantData The bgen data on this variant.
     * @param phenoName The phenotype name.
     * @param phenotypes The phenotype values.
     * @param nValidValues The number of phenotypes that are finite and not NA.
     */
    private void runLinearModel(
            int variantIndex,
            int[] testedAlleleIndexes,
            BgenVariantData variantData,
            String phenoName,
            double[] phenotypes
    ) {

        VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

        for (int alleleI : testedAlleleIndexes) {

            int[] childIndexes = covariatesHandler.originalIndexMap.get(phenoName);

            // Check transmission
            double cMin = Double.NaN;
            double cMax = Double.NaN;
            double mMin = Double.NaN;
            double mMax = Double.NaN;
            double fMin = Double.NaN;
            double fMax = Double.NaN;

            for (int childIndex : childIndexes) {

                String childId = childToParentMap.children[childIndex];
                String motherId = childToParentMap.getMother(childId);
                String fatherId = childToParentMap.getFather(childId);

                if (variantData.contains(childId)) {

                    double cValue = variantData.getSummedProbability(childId, alleleI);

                    if (Double.isNaN(cMin) || cValue < cMin) {

                        cMin = cValue;

                    }
                    if (Double.isNaN(cMax) || cValue > cMax) {

                        cMax = cValue;

                    }
                }

                if (variantData.contains(motherId)) {

                    double mValue = variantData.getSummedProbability(motherId, alleleI);

                    if (Double.isNaN(mMin) || mValue < mMin) {

                        mMin = mValue;

                    }
                    if (Double.isNaN(mMax) || mValue > mMax) {

                        mMax = mValue;

                    }
                }

                if (variantData.contains(fatherId)) {

                    double fValue = variantData.getSummedProbability(fatherId, alleleI);

                    if (Double.isNaN(fMin) || fValue < fMin) {

                        fMin = fValue;

                    }
                    if (Double.isNaN(fMax) || fValue > fMax) {

                        fMax = fValue;

                    }
                }
            }

            if (!x0 || cMax - cMin > 0.5 || mMax - mMin > 0.5 || fMax - fMin > 0.5) {

                // Prepare the objects to use for the models, gather values for the histograms
                double[] phenoValues = phenotypesHandler.phenoMap.get(phenoName);
                double phenoMean = phenotypesHandler.phenoMeanMap.get(phenoName);
                double[] rss0s = new double[phenoValues.length];

                ArrayList<RegressionResult> regressionResults = new ArrayList<>(models.length);

                for (int i = 0; i < models.length; i++) {

                    Model model = models[i];

                    regressionResults.add(
                            new RegressionResult(
                                    model,
                                    phenoValues.length
                            )
                    );
                }

                TreeMap<Integer, Integer>[] hHist = new TreeMap[4];

                for (int j = 0; j < 4; j++) {

                    hHist[j] = new TreeMap<>();

                }

                TreeMap<Integer, Integer> childHist = new TreeMap<>();
                TreeMap<Integer, Integer> motherHist = new TreeMap<>();
                TreeMap<Integer, Integer> fatherHist = new TreeMap<>();

                for (int i = 0; i < childIndexes.length; i++) {

                    int childIndex = childIndexes[i];

                    String childId = childToParentMap.children[childIndex];
                    String motherId = childToParentMap.getMother(childId);
                    String fatherId = childToParentMap.getFather(childId);

                    if (variantData.contains(childId) && variantData.contains(motherId) && variantData.contains(fatherId)) {

                        double[] haplotypes = variantData.getHaplotypes(
                                childId,
                                motherId,
                                fatherId,
                                alleleI
                        );

                        double distY = phenotypes[i] - phenoMean;
                        rss0s[i] = distY * distY;

                        for (int j = 0; j < 4; j++) {

                            int hJ = (int) Math.round(haplotypes[j]);

                            TreeMap<Integer, Integer> hHistJ = hHist[j];
                            Integer frequency = hHistJ.get(hJ);

                            if (frequency != null) {

                                hHistJ.put(hJ, frequency + 1);

                            } else {

                                hHistJ.put(hJ, 1);

                            }
                        }
                    }

                    if (variantData.contains(childId)) {

                        int nAltChild = (int) Math.round(variantData.getSummedProbability(childId, alleleI));

                        Integer frequency = childHist.get(nAltChild);

                        if (frequency != null) {

                            childHist.put(nAltChild, frequency + 1);

                        } else {

                            childHist.put(nAltChild, 1);

                        }
                    }

                    if (variantData.contains(motherId)) {

                        int nAltMother = (int) Math.round(variantData.getSummedProbability(motherId, alleleI));

                        Integer frequency = motherHist.get(nAltMother);

                        if (frequency != null) {

                            motherHist.put(nAltMother, frequency + 1);

                        } else {

                            motherHist.put(nAltMother, 1);

                        }
                    }

                    if (variantData.contains(fatherId)) {

                        int nAltFather = (int) Math.round(variantData.getSummedProbability(fatherId, alleleI));

                        Integer frequency = fatherHist.get(nAltFather);

                        if (frequency != null) {

                            fatherHist.put(nAltFather, frequency + 1);

                        } else {

                            fatherHist.put(nAltFather, 1);

                        }
                    }
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
                        "hmnt(",
                        getHistogramAsString(hHist[0]),
                        ");hmt(",
                        getHistogramAsString(hHist[1]),
                        ");hft(",
                        getHistogramAsString(hHist[2]),
                        ");hfnt(",
                        getHistogramAsString(hHist[3]),
                        ")"
                );

                // Get matrices for haplotypes and individuals
                double[][] haplotypeX = new double[childIndexes.length][4];
                double[][] childX = new double[childIndexes.length][1];
                double[][] motherX = new double[childIndexes.length][1];
                double[][] fatherX = new double[childIndexes.length][1];

                for (int i = 0; i < childIndexes.length; i++) {

                    int childIndex = childIndexes[i];

                    String childId = childToParentMap.children[childIndex];
                    String motherId = childToParentMap.getMother(childId);
                    String fatherId = childToParentMap.getFather(childId);

                    if (variantData.contains(childId)) {

                        double[] haplotypes = variantData.getHaplotypes(
                                childId,
                                motherId,
                                fatherId,
                                alleleI
                        );
                        haplotypeX[i][0] = haplotypes[0];
                        haplotypeX[i][1] = haplotypes[1];
                        haplotypeX[i][2] = haplotypes[2];
                        haplotypeX[i][3] = haplotypes[3];

                        childX[i][0] = variantData.getSummedProbability(childId, alleleI);

                    }

                    if (variantData.contains(motherId)) {

                        motherX[i][0] = variantData.getSummedProbability(motherId, alleleI);

                    }

                    if (variantData.contains(fatherId)) {

                        fatherX[i][0] = variantData.getSummedProbability(fatherId, alleleI);

                    }
                }

                // Adjust for covariates
                haplotypeX = covariatesHandler.getAdjustedValues(phenoName, haplotypeX);
                childX = covariatesHandler.getAdjustedValues(phenoName, childX);
                motherX = covariatesHandler.getAdjustedValues(phenoName, motherX);
                fatherX = covariatesHandler.getAdjustedValues(phenoName, fatherX);

                // Run the regressions
                OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
                HashMap<String, RegressionResult> regressionRestultsMap = new HashMap<>(regressionResults.size());

                for (int modelI = 0; modelI < models.length; modelI++) {

                    Model model = models[modelI];

                    // Get the samples to use
                    ArrayList<Integer> rowsToRun = new ArrayList<>(childIndexes.length);

                    for (int i = 0; i < childIndexes.length; i++) {

                        int childIndex = childIndexes[i];

                        String childId = childToParentMap.children[childIndex];
                        String motherId = childToParentMap.getMother(childId);
                        String fatherId = childToParentMap.getFather(childId);

                        if (Model.hasData(model, childId, motherId, fatherId, variantData)) {

                            rowsToRun.add(i);

                        }
                    }

                    // Create the matrices to use in the regression.
                    double[][] x = new double[rowsToRun.size()][model.betaNames.length];
                    double[] min = new double[model.betaNames.length];
                    Arrays.fill(min, Double.NaN);
                    double[] max = new double[model.betaNames.length];
                    Arrays.fill(max, Double.NaN);
                    double[] y = new double[rowsToRun.size()];
                    double rss0 = 0.0;

                    for (int i = 0; i < rowsToRun.size(); i++) {

                        int childIndex = rowsToRun.get(i);

                        y[i] = phenoValues[childIndex];
                        rss0 += rss0s[childIndex];

                        for (int j = 0; j < model.betaNames.length; j++) {

                            double xValue = Model.getXValueAt(
                                    model,
                                    childIndex,
                                    j,
                                    haplotypeX,
                                    childX,
                                    motherX,
                                    fatherX
                            );

                            x[i][j] = xValue;

                            double minValue = min[j];

                            if (Double.isNaN(minValue) || xValue < minValue) {

                                min[j] = xValue;

                            }

                            double maxValue = max[j];

                            if (Double.isNaN(maxValue) || xValue > maxValue) {

                                max[j] = xValue;

                            }
                        }
                    }

                    // Check singularities.
                    boolean singularity = false;

                    for (int betaI = 0; betaI < min.length; betaI++) {

                        if (Double.isNaN(max[betaI]) || Double.isNaN(min[betaI]) || max[betaI] - min[betaI] < 0.1) {

                            singularity = true;

                        }
                    }

                    if (!singularity) {

                        // Run regression
                        RegressionResult regressionResult = regressionResults.get(modelI);

                        try {

                            regression.newSampleData(y, x);
                            double[] betas = regression.estimateRegressionParameters();
                            double[] betaStandardErrors = regression.estimateRegressionParametersStandardErrors();
                            double[] betaResiduals = regression.estimateResiduals();

                            regressionResult.beta = Arrays.copyOfRange(betas, 1, betas.length);
                            regressionResult.betaStandardError = Arrays.copyOfRange(betaStandardErrors, 1, betaStandardErrors.length);
                            regressionResult.betaResiduals = Arrays.copyOfRange(betaResiduals, 1, betaResiduals.length);

                            regressionResult.computeRSS();
                            regressionResult.computeModelSignificance(rss0);
                            regressionResult.computeBetaSignificance();

                            regressionRestultsMap.put(model.name(), regressionResult);

                        } catch (SingularMatrixException singularMatrixException) {

                            logger.logVariant(
                                    variantInformation.id,
                                    String.join(" ",
                                            model.name(),
                                            phenoName,
                                            "Singularity detected",
                                            "SinglularityException"
                                    )
                            );

                            if (debugSingularities) {

                                writeSingularityDebugReport(
                                        variantInformation.id,
                                        phenoName,
                                        model.name(),
                                        x,
                                        phenoValues
                                );
                            }
                        }
                    } else {

                        logger.logVariant(
                                variantInformation.id,
                                String.join(" ",
                                        model.name(),
                                        phenoName,
                                        "Singularity anticipated",
                                        "xmax - xmin < 0.5"
                                )
                        );
                    }
                }

                // Estimate model significance relative to parent models
                regressionRestultsMap.values()
                        .forEach(
                                regressionResult -> regressionResult.computeModelSignificance(
                                        regressionRestultsMap
                                )
                        );

                // Estimate the share of mendelian errors
                double mendelianErrors = MendelianErrorEstimator.estimateMendelianErrorPrevalence(variantData, childToParentMap, alleleI);

                // Export
                StringBuilder stringBuilder = new StringBuilder();
                stringBuilder
                        .append(phenoName)
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.contig)
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.position)
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.id)
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.rsid)
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.alleles[alleleI])
                        .append(IoUtils.SEPARATOR)
                        .append(variantInformation.getOtherAllele(alleleI))
                        .append(IoUtils.SEPARATOR)
                        .append(childIndexes.length)
                        .append(IoUtils.SEPARATOR)
                        .append(altHistograms)
                        .append(IoUtils.SEPARATOR)
                        .append(hHistograms)
                        .append(IoUtils.SEPARATOR)
                        .append(mendelianErrors);

                regressionResults
                        .forEach(
                                regressionResult -> regressionResult.appendResults(stringBuilder)
                        );

                String line = stringBuilder
                        .append(IoUtils.LINE_SEPARATOR)
                        .toString();

                gzIndexMutex.acquire();

                IndexedGzCoordinates coordinates = outputWriter.append(line);

                resultsIndex.writeLine(variantInformation.contig,
                        Integer.toString(variantInformation.position),
                        variantInformation.id,
                        variantInformation.rsid,
                        phenoName,
                        Integer.toString(coordinates.compressedLength),
                        Integer.toString(coordinates.uncompressedLength)
                );

                gzIndexMutex.release();

            } else {

                logger.logVariant(
                        variantInformation.id,
                        "Same alleles in all individuals."
                );
            }
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
                                Double.toString(entry.getKey()),
                                Integer.toString(entry.getValue())
                        )
                )
                .collect(
                        Collectors.joining(",")
                );
    }

    /**
     * Writes a debug report when a singularity is found.
     *
     * @param variantId The id of the variant.
     * @param phenoName The name of the phenotype.
     * @param modelName The name of the model.
     * @param x The x values.
     * @param y The y values.
     */
    private void writeSingularityDebugReport(
            String variantId,
            String phenoName,
            String modelName,
            double[][] x,
            double[] y
    ) {

        String debugFileName = String.join("_", variantId, phenoName, modelName, "debug_singularity_x");
        File singularityDebugFile = new File(outputWriter.getFile().getParentFile(), debugFileName);

        try (SimpleFileWriter writer = new SimpleFileWriter(singularityDebugFile, true)) {

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

        try (SimpleFileWriter writer = new SimpleFileWriter(singularityDebugFile, true)) {

            for (int yi = 0; yi < y.length; yi++) {

                writer.writeLine(Double.toString(y[yi]));

            }
        }
    }
}
