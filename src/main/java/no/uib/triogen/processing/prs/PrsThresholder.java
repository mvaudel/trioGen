package no.uib.triogen.processing.prs;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.utils.SimpleSemaphore;
import no.uib.triogen.utils.Utils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 * This class scores genotypes based on pruned results from PrsTrainer.
 *
 * @author Marc Vaudel
 */
public class PrsThresholder {

    /**
     * The lowest p-value threshold.
     */
    public static final double LOWEST_P_VALUE = 1e-6;
    /**
     * The step to use for the p-value threshold decrease. 2 means divided by 2.
     */
    public static final double P_VALUE_STEP = 2;
    /**
     * The beta quantiles to try.
     */
    public static final double[] BETA_QUANTILES = new double[]{0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    /**
     * The folder to use to store the bgen index files.
     */
    private final File bgenIndexFolder;
    /**
     * The path top the file containing the genotypes.
     */
    private final String genotypesFilePath;
    /**
     * The name of the phenotype.
     */
    private final String phenoName;
    /**
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The file to export the results to.
     */
    private final File destinationFile;
    /**
     * The file containing the training summary stats.
     */
    private final File trainingFile;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The pattern to use to find the effect size column for each variable in
     * the model.
     */
    private final String betaColumnPattern;
    /**
     * The pattern to use to find the standard error column for each variable in
     * the model.
     */
    private final String seColumnPattern;
    /**
     * The trio model to use.
     */
    private final Model model;
    /**
     * The ordered names of the variables used in the model.
     */
    private final String[] variableNames;
    /**
     * The number of bins to use.
     */
    private final int nBins;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * Semaphore to synchronize the scoring threads.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param genotypesFilePath The path to the file containing the genotypes.
     * @param phenotypesFile The file containing the phenotypes.
     * @param phenoName The name of the phenotype to process.
     * @param childToParentMap The map of trios.
     * @param trainingFile The file containing the training data.
     * @param destinationFile The file to export the result to.
     * @param variantList The targets to use.
     * @param model The trio model to use.
     * @param variableNames The names of the variables to include.
     * @param betaColumnPattern The pattern to use to find the effect size
     * column for each variable in the model.
     * @param seColumnPattern The pattern to use to find the standard error
     * column for each variable in the model.
     * @param nBins The number of bins to use.
     * @param bgenIndexFolder The folder for the bgen index files.
     * @param logger The logger.
     */
    public PrsThresholder(
            String genotypesFilePath,
            File phenotypesFile,
            String phenoName,
            ChildToParentMap childToParentMap,
            File trainingFile,
            File destinationFile,
            VariantList variantList,
            String betaColumnPattern,
            String seColumnPattern,
            Model model,
            String[] variableNames,
            int nBins,
            File bgenIndexFolder,
            SimpleCliLogger logger
    ) {

        this.genotypesFilePath = genotypesFilePath;
        this.phenotypesFile = phenotypesFile;
        this.phenoName = phenoName;
        this.childToParentMap = childToParentMap;
        this.trainingFile = trainingFile;
        this.destinationFile = destinationFile;
        this.variantList = variantList;
        this.betaColumnPattern = betaColumnPattern;
        this.seColumnPattern = seColumnPattern;
        this.model = model;
        this.variableNames = variableNames;
        this.nBins = nBins;
        this.bgenIndexFolder = bgenIndexFolder;
        this.logger = logger;

    }

    /**
     * Runs the scoring.
     *
     * @throws IOException Exception thrown if an error occurs while reading or
     * writing a file.
     */
    public void run() throws IOException {

        logger.logMessage("Thresholding risk score");

        long allStart = Instant.now().getEpochSecond();

        logger.logMessage("Loading phenotypes");

        long start = Instant.now().getEpochSecond();

        HashSet<String> phenoNames = new HashSet<>(1);
        phenoNames.add(phenoName);

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                phenoNames
        );

        double[] phenotypes = phenotypesHandler.phenoMap.get(phenoName);

        double[] binnedPhenotypes = Utils.bin(phenotypes, nBins);

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Loading phenotypes done (" + duration + " seconds).");

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            writer.writeLine(
                    "scoring_mode",
                    "pValue_threshold",
                    "beta_quantile",
                    "bin_r2"
            );

            for (PrsUtils.ScoringMode scoringMode : PrsUtils.ScoringMode.values()) {

                double pValueThreshold = 0.05;

                while (pValueThreshold >= LOWEST_P_VALUE) {

                    HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData = PrsUtils.parseScoringData(
                            trainingFile,
                            betaColumnPattern,
                            seColumnPattern,
                            variableNames,
                            variantList,
                            pValueThreshold,
                            scoringMode
                    );

                    for (double betaQuantile : BETA_QUANTILES) {

                        logger.logMessage("Thresholding {scoringMode: " + scoringMode + ", pValueThreshold: " + pValueThreshold + ", betaQuantile: " + betaQuantile + "}");

                        start = Instant.now().getEpochSecond();

                        NormalDistribution normalDistribution = new NormalDistribution(0, 1);
                        double seScaling = normalDistribution.inverseCumulativeProbability(betaQuantile);

                        HashMap<String, double[]> scores = PrsComputer.getScores(
                                genotypesFilePath,
                                scoringData,
                                childToParentMap,
                                model,
                                variableNames,
                                seScaling,
                                bgenIndexFolder,
                                semaphore,
                                logger
                        );

                        double[] scoresSum = new double[childToParentMap.children.length];

                        for (int i = 0; i < childToParentMap.children.length; i++) {

                            String childId = childToParentMap.children[i];

                            double[] childScores = scores.get(childId);

                            double childScore = Arrays.stream(childScores).sum();

                            scoresSum[i] = childScore;

                        }

                        double[] binnedScores = Utils.bin(scoresSum, nBins);

                        PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();

                        double correlation = pearsonsCorrelation.correlation(binnedScores, binnedPhenotypes);

                        writer.writeLine(
                                scoringMode.name(),
                                Double.toString(pValueThreshold),
                                Double.toString(betaQuantile),
                                Double.toString(correlation)
                        );

                    }
                    
                    pValueThreshold /= 2;
                    
                }
            }
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Scoring done and exported to " + destinationFile + " (" + duration + " seconds).");

        long allEnd = Instant.now().getEpochSecond();
        long allDuration = allEnd - allStart;

        logger.logMessage("Computing risk score done (" + allDuration + " seconds).");

    }
}
