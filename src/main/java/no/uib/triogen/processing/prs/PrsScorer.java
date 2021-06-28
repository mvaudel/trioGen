package no.uib.triogen.processing.prs;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.HashMap;
import java.util.stream.IntStream;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.trio_genotypes.VariantList;
import static no.uib.triogen.processing.prs.PrsTrainer.VARIABLE_WILDCARD;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class scores genotypes based on pruned results from PrsTrainer.
 *
 * @author Marc Vaudel
 */
public class PrsScorer {

    /**
     * Wildcard for the chromosome name in the genotypes file.
     */
    public static final String CHROMOSOME_WILDCARD = "{chr}";
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
     * @param pValueThreshold The p-value threshold to use.
     * @param scorerMode The scorer mode to use.
     * @param logger The logger.
     */
    public PrsScorer(
            String genotypesFilePath,
            ChildToParentMap childToParentMap,
            File trainingFile,
            File destinationFile,
            VariantList variantList,
            String betaColumnPattern,
            Model model,
            String[] variableNames,
            double pValueThreshold,
            ScoringMode scorerMode,
            SimpleCliLogger logger
    ) {

        this.genotypesFilePath = genotypesFilePath;
        this.childToParentMap = childToParentMap;
        this.trainingFile = trainingFile;
        this.destinationFile = destinationFile;
        this.variantList = variantList;
        this.betaColumnPattern = betaColumnPattern;
        this.model = model;
        this.variableNames = variableNames;
        this.pValueThreshold = pValueThreshold;
        this.scoringMode = scorerMode;
        this.logger = logger;

    }

    /**
     * Runs the scoring.
     *
     * @throws IOException Exception thrown if an error occurs while reading or
     * writing a file.
     */
    public void run() throws IOException {

        logger.logMessage("Parsing scoring data from " + trainingFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        HashMap<String, HashMap<String, HashMap<String, double[]>>> scoringData = parseScoringData();

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing scoring data from " + trainingFile + " done (" + duration + " seconds), " + scoringData.size() + " variants used for scoring.");

        logger.logMessage("Scoring");

        start = Instant.now().getEpochSecond();

        HashMap<String, double[]> scores = getScores(scoringData);

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

                                writer.writeLine(header.toString());

                            }
                    );
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Export done (" + duration + " seconds).");

    }

    private HashMap<String, double[]> getScores(HashMap<String, HashMap<String, HashMap<String, double[]>>> scoringData) {

        HashMap<String, double[]> scores = new HashMap<>(childToParentMap.children.length);

        for (String childId : childToParentMap.children) {

            scores.put(childId, new double[variableNames.length]);

        }

        scoringData.entrySet()
                .stream()
                .parallel()
                .forEach(
                        entry -> processChromosome(
                                entry.getKey(),
                                entry.getValue(),
                                scores
                        )
                );

        return scores;

    }

    private void processChromosome(
            String chromosome,
            HashMap<String, HashMap<String, double[]>> chromosomeScoringData,
            HashMap<String, double[]> scores
    ) {

        try {

            File genotypesFile = new File(genotypesFilePath.replace(CHROMOSOME_WILDCARD, chromosome));

            logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

            long start = Instant.now().getEpochSecond();

        HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(chromosome);

        if (inheritanceMap == null) {

            throw new IllegalArgumentException("Mode of inheritance not implemented for " + chromosome + ".");

        }

        int defaultMotherPloidy = InheritanceUtils.getDefaultMotherPloidy(chromosome);
        int defaultFatherPloidy = InheritanceUtils.getDefaultFatherPloidy(chromosome);

            BgenIndex bgenIndex = BgenIndex.getBgenIndex(genotypesFile);
            BgenFileReader bgenFileReader = new BgenFileReader(
                    genotypesFile,
                    bgenIndex,
                    inheritanceMap,
                    defaultMotherPloidy,
                    defaultFatherPloidy
            );

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            logger.logMessage("Parsing " + genotypesFile + " done (" + duration + " seconds)");

            logger.logMessage("Scoring chromosome " + chromosome);

            start = Instant.now().getEpochSecond();

            IntStream.range(0, bgenIndex.variantIdArray.length)
                    .parallel()
                    .filter(
                            variantIndex -> chromosomeScoringData.containsKey(bgenIndex.variantIdArray[variantIndex])
                    )
                    .forEach(
                            variantIndex -> score(
                                    variantIndex,
                                    chromosomeScoringData.get(bgenIndex.variantIdArray[variantIndex]),
                                    bgenFileReader,
                                    bgenIndex,
                                    scores
                            )
                    );

            end = Instant.now().getEpochSecond();
            duration = end - start;

            logger.logMessage("Scoring chromosome " + chromosome + " done (" + duration + " seconds)");

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    private void score(
            int variantIndex,
            HashMap<String, double[]> variantScoringData,
            BgenFileReader bgenFileReader,
            BgenIndex bgenIndex,
            HashMap<String, double[]> scores
    ) {

        variantScoringData.entrySet()
                .forEach(
                        entry -> score(
                                variantIndex,
                                entry.getKey(),
                                entry.getValue(),
                                bgenFileReader,
                                bgenIndex,
                                scores
                        )
                );

    }

    private void score(
            int variantIndex,
            String effectAllele,
            double[] scoringDetails,
            BgenFileReader bgenFileReader,
            BgenIndex bgenIndex,
            HashMap<String, double[]> scores
    ) {

        // Scoring details
        double weight = scoringDetails[0];

        // Get the index of the allele to score
        VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

        int alleleI = -1;

        for (int i = 0; i < variantInformation.alleles.length; i++) {

            if (variantInformation.alleles[i].equals(effectAllele)) {

                alleleI = i;
                break;

            }
        }

        if (alleleI == -1) {

            throw new IllegalArgumentException("Allele " + effectAllele + " not found for variant " + variantInformation.id);

        }

        // Parse genotypes

            ZstdDecompressor decompressor = new ZstdDecompressor();
            
        BgenVariantData variantData = bgenFileReader.getVariantData(variantIndex);
        variantData.parse(
                childToParentMap,
                decompressor
        );

        // Get matrices for haplotypes and individuals
        double[][] haplotypeX = new double[childToParentMap.children.length][4];
        double[][] childX = new double[childToParentMap.children.length][1];
        double[][] motherX = new double[childToParentMap.children.length][1];
        double[][] fatherX = new double[childToParentMap.children.length][1];

        for (int childI = 0; childI < childToParentMap.children.length; childI++) {

            String childId = childToParentMap.children[childI];
            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            if (variantData.contains(childId)) {

                double[] haplotypes = variantData.getHaplotypes(
                        childId,
                        motherId,
                        fatherId,
                        alleleI
                );
                haplotypeX[childI][0] = haplotypes[0];
                haplotypeX[childI][1] = haplotypes[1];
                haplotypeX[childI][2] = haplotypes[2];
                haplotypeX[childI][3] = haplotypes[3];

                childX[childI][0] = variantData.getSummedProbability(childId, alleleI);

            }

            if (variantData.contains(motherId)) {

                motherX[childI][0] = variantData.getSummedProbability(motherId, alleleI);

            }

            if (variantData.contains(fatherId)) {

                fatherX[childI][0] = variantData.getSummedProbability(fatherId, alleleI);

            }

            semaphore.acquire();

            double[] childScores = scores.get(childId);

            for (int variableI = 0; variableI < model.betaNames.length; variableI++) {

                double xValue = Model.getXValueAt(
                        model,
                        childI,
                        variableI,
                        haplotypeX,
                        childX,
                        motherX,
                        fatherX
                );

                double beta = scoringDetails[variableI + 1];

                double variantContribution;

                switch (scoringMode) {

                    case lead:

                        variantContribution = xValue * beta;
                        break;

                    case weighted:
                        variantContribution = xValue * beta / weight;
                        break;

                    default:
                        throw new UnsupportedOperationException("Scoring mode " + scoringMode + " not implemented.");

                }

                childScores[variableI] = childScores[variableI] + variantContribution;

            }

            semaphore.release();

        }
    }

    /**
     * For each variant, parses the weight and effect size.
     *
     * @return The scoring data in a map, chromosome to variant id to effect
     * allele to weight and effect size.
     */
    private HashMap<String, HashMap<String, HashMap<String, double[]>>> parseScoringData() {

        HashMap<String, HashMap<String, HashMap<String, double[]>>> scoringData = new HashMap<>();

        int[] betaColumnIndexes = new int[variableNames.length];

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(trainingFile)) {

            String line = reader.readLine();

            String[] lineSplit = line.split(SEPARATOR);

            for (int i = 0; i < lineSplit.length; i++) {

                for (int j = 0; j < variableNames.length; j++) {

                    String variable = variableNames[j];

                    String betaColumn = betaColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(betaColumn)) {

                        betaColumnIndexes[j] = i;

                    }
                }
            }

            for (int j = 0; j < variableNames.length; j++) {

                String variable = variableNames[j];

                if (betaColumnIndexes[j] == -1) {

                    throw new IllegalArgumentException("Effect size column not found for '" + variable + "'.");

                }
            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(SEPARATOR);

                String leadVariantId = lineSplit[0];
                double leadP = Double.parseDouble(lineSplit[2]);

                String variantId = lineSplit[3];
                String variantRsid = lineSplit[4];

                if (variantList != null && variantList.contains(variantId)
                        || variantList != null && variantList.contains(variantRsid) 
                        || variantList == null && leadP <= pValueThreshold) {

                    if (scoringMode == ScoringMode.weighted || leadVariantId.equals(variantId)) {

                        String chromosome = lineSplit[5];
                        String ea = lineSplit[8];

                        double weight = Double.parseDouble(lineSplit[9]);

                        double[] result = new double[variableNames.length + 1];

                        result[0] = weight;

                        for (int j = 0; j < variableNames.length; j++) {

                            String betaString = lineSplit[betaColumnIndexes[j]];

                            double beta = Double.parseDouble(betaString);

                            result[j + 1] = beta;

                        }

                        HashMap<String, HashMap<String, double[]>> chromosomeValues = scoringData.get(chromosome);

                        if (chromosomeValues == null) {

                            chromosomeValues = new HashMap<>();
                            scoringData.put(chromosome, chromosomeValues);

                        }

                        HashMap<String, double[]> variantValues = chromosomeValues.get(variantId);

                        if (variantValues == null) {

                            variantValues = new HashMap<>(1);
                            chromosomeValues.put(variantId, variantValues);

                        }

                        variantValues.put(ea, result);

                    }
                }
            }
        }

        return scoringData;

    }

    public enum ScoringMode {

        lead(0),
        weighted(1);

        private ScoringMode(int index) {

            this.index = index;

        }

        public final int index;

    }
}
