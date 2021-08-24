package no.uib.triogen.processing.prs;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.processing.prs.PrsUtils.ScoringMode;
import no.uib.triogen.utils.SimpleSemaphore;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * This class scores genotypes based on pruned results from PrsTrainer.
 *
 * @author Marc Vaudel
 */
public class PrsScorer {

    /**
     * The folder to use to store the bgen index files.
     */
    private final File bgenIndexFolder;
    /**
     * The path top the file containing the genotypes.
     */
    private final String genotypesFilePath;
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
     * The lead p-value threshold.
     */
    private final double pValueThreshold;
    /**
     * The scaling value for the standard error.
     */
    private final double seScaling;
    /**
     * The scorer mode to use.
     */
    private final ScoringMode scoringMode;
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
     * @param pValueThreshold The p-value threshold to use.
     * @param qBeta The quantile to use on the beta estimation.
     * @param scorerMode The scorer mode to use.
     * @param bgenIndexFolder The folder for the bgen index files.
     * @param logger The logger.
     */
    public PrsScorer(
            String genotypesFilePath,
            ChildToParentMap childToParentMap,
            File trainingFile,
            File destinationFile,
            VariantList variantList,
            String betaColumnPattern,
            String seColumnPattern,
            Model model,
            String[] variableNames,
            double pValueThreshold,
            double qBeta,
            ScoringMode scorerMode,
            File bgenIndexFolder,
            SimpleCliLogger logger
    ) {

        this.genotypesFilePath = genotypesFilePath;
        this.childToParentMap = childToParentMap;
        this.trainingFile = trainingFile;
        this.destinationFile = destinationFile;
        this.variantList = variantList;
        this.betaColumnPattern = betaColumnPattern;
        this.seColumnPattern = seColumnPattern;
        this.model = model;
        this.variableNames = variableNames;
        this.pValueThreshold = pValueThreshold;
        this.scoringMode = scorerMode;
        this.bgenIndexFolder = bgenIndexFolder;
        this.logger = logger;

        NormalDistribution normalDistribution = new NormalDistribution(0, 1);
        this.seScaling = normalDistribution.inverseCumulativeProbability(qBeta);

    }

    /**
     * Runs the scoring.
     *
     * @throws IOException Exception thrown if an error occurs while reading or
     * writing a file.
     */
    public void run() throws IOException {

        logger.logMessage("Computing risk score");

        long allStart = Instant.now().getEpochSecond();

        logger.logMessage("Parsing scoring data from " + trainingFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData = PrsUtils.parseScoringData(
                trainingFile, 
                betaColumnPattern, 
                seColumnPattern, 
                variableNames, 
                variantList, 
                pValueThreshold, 
                scoringMode
        );

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing scoring data from " + trainingFile + " done (" + duration + " seconds).");

        logger.logMessage("Scoring");

        start = Instant.now().getEpochSecond();

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

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Scoring done (" + duration + " seconds).");

        logger.logMessage("Exporting to " + destinationFile);

        start = Instant.now().getEpochSecond();

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            StringBuilder header = new StringBuilder();
            header.append("child_id");

            for (String variable : variableNames) {

                header.append(IoUtils.SEPARATOR)
                        .append(variable);

            }

            writer.writeLine(header.toString());

            scores.entrySet()
                    .stream()
                    .forEach(
                            entry -> {

                                StringBuilder line = new StringBuilder();
                                line.append(entry.getKey());

                                for (int variableI = 0; variableI < variableNames.length; variableI++) {

                                    line.append(IoUtils.SEPARATOR)
                                            .append(Double.toString(entry.getValue()[variableI]));

                                }

                                writer.writeLine(line.toString());

                            }
                    );
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Export done (" + duration + " seconds).");

        long allEnd = Instant.now().getEpochSecond();
        long allDuration = allEnd - allStart;

        logger.logMessage("Computing risk score done (" + allDuration + " seconds).");

    }
}
