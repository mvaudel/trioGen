package no.uib.triogen.processing.prs;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
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
import no.uib.triogen.utils.Utils;
import org.apache.commons.math3.distribution.NormalDistribution;

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

        logger.logMessage("Parsing scoring data from " + trainingFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData = parseScoringData();

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing scoring data from " + trainingFile + " done (" + duration + " seconds).");

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

                                writer.writeLine(line.toString());

                            }
                    );
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Export done (" + duration + " seconds).");

    }

    private HashMap<String, double[]> getScores(HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData) {

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
            HashMap<String, HashMap<String, HashMap<String, double[]>>> chromosomeScoringData,
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

            logger.logMessage("Gathering genotyped variants from " + chromosome);

            int nVariantsScore = 0;
            int nLociScore = chromosomeScoringData.size();

            start = Instant.now().getEpochSecond();

            HashMap<String, String> scoringToLeadVariantMap = new HashMap<>();

            for (Entry<String, HashMap<String, HashMap<String, double[]>>> entry : chromosomeScoringData.entrySet()) {

                String leadVariant = entry.getKey();

                for (String scoringVariant : entry.getValue().keySet()) {

                    if (scoringToLeadVariantMap.containsKey(scoringVariant)) {

                        String otherleadVariant = scoringToLeadVariantMap.get(scoringVariant);

                        throw new IllegalArgumentException(scoringVariant + " found in two lead variants: " + leadVariant + " and " + otherleadVariant + ".");

                    }

                    scoringToLeadVariantMap.put(scoringVariant, leadVariant);

                    nVariantsScore++;

                }
            }

            HashMap<String, HashMap<String, Integer>> leadToScoringVariantIndexMap = new HashMap<>();
            int nVariants = 0;

            for (int bgenI = 0; bgenI < bgenIndex.variantIdArray.length; bgenI++) {

                String bgenVariantId = bgenIndex.variantIdArray[bgenI];

                String leadVariant = scoringToLeadVariantMap.get(bgenVariantId);

                if (leadVariant != null) {

                    HashMap<String, Integer> scoringVariants = leadToScoringVariantIndexMap.get(leadVariant);

                    if (scoringVariants == null) {

                        scoringVariants = new HashMap<>(1);
                        leadToScoringVariantIndexMap.put(leadVariant, scoringVariants);

                    }

                    scoringVariants.put(bgenVariantId, bgenI);

                    nVariants++;

                }
            }

            int nLociFound = leadToScoringVariantIndexMap.size();

            end = Instant.now().getEpochSecond();
            duration = end - start;

            logger.logMessage("Gathering genotyped variants from " + chromosome + " done (" + duration + " seconds), " + nVariants + "/" + nLociFound + " variants/loci found of " + nVariantsScore + "/" + nLociScore + ".");

            logger.logMessage("Scoring chromosome " + chromosome);

            start = Instant.now().getEpochSecond();

            scoringToLeadVariantMap.entrySet()
                    .stream()
                    .parallel()
                    .forEach(
                            entry -> {

                                String scoringVariant = entry.getKey();
                                String leadVariant = entry.getValue();
                                HashMap<String, Integer> indexMap = leadToScoringVariantIndexMap.get(leadVariant);
                                int bgenVariantIndex = indexMap.get(scoringVariant);
                                double weight = 1.0 / ((double) indexMap.size());
                                HashMap<String, double[]> alleles = chromosomeScoringData.get(leadVariant).get(scoringVariant);

                                for (Entry<String, double[]> allelesEntry : alleles.entrySet()) {

                                    String allele = allelesEntry.getKey();
                                    double[] scoringDetails = allelesEntry.getValue();

                                    score(bgenVariantIndex, allele, weight, scoringDetails, bgenFileReader, bgenIndex, scores);

                                }
                            }
                    );

            end = Instant.now().getEpochSecond();
            duration = end - start;

            logger.logMessage("Scoring chromosome " + chromosome + " done (" + duration + " seconds)");

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    private void score(
            int bgenVariantIndex,
            String effectAllele,
            double weight,
            double[] scoringDetails,
            BgenFileReader bgenFileReader,
            BgenIndex bgenIndex,
            HashMap<String, double[]> scores
    ) {

        // Get the index of the allele to score
        VariantInformation variantInformation = bgenIndex.variantInformationArray[bgenVariantIndex];

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

        BgenVariantData variantData = bgenFileReader.getVariantData(bgenVariantIndex);
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

                double beta = scoringDetails[variableI];
                double se = scoringDetails[variableI + variableNames.length];

                double betaEstimate = 0.0;

                if (beta > 0.0) {

                    betaEstimate = beta + seScaling * se;

                    if (betaEstimate < 0.0) {

                        betaEstimate = 0.0;

                    }
                } else if (beta < 0.0) {

                    betaEstimate = beta - seScaling * se;

                    if (betaEstimate > 0.0) {

                        betaEstimate = 0.0;

                    }
                }

                double variantContribution = xValue * betaEstimate * weight;

                childScores[variableI] = childScores[variableI] + variantContribution;

            }

            semaphore.release();

        }
    }

    /**
     * For each variant, parses the weight and effect size.
     *
     * @return The scoring data in a map, chromosome to lead variant to scoring
     * variant to effect allele to weight and effect size.
     */
    private HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> parseScoringData() {

        HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData = new HashMap<>();

        int[] betaColumnIndexes = new int[variableNames.length];
        Arrays.fill(betaColumnIndexes, -1);
        
        int[] seColumnIndexes = new int[variableNames.length];
        Arrays.fill(seColumnIndexes, -1);

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

                    String seColumn = seColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(seColumn)) {

                        seColumnIndexes[j] = i;

                    }
                }
            }

            for (int j = 0; j < variableNames.length; j++) {

                String variable = variableNames[j];

                if (betaColumnIndexes[j] == -1) {

                    String betaColumn = betaColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);
                    throw new IllegalArgumentException("Effect size column '" + betaColumn + "' not found for '" + variable + "'.");

                }

                if (seColumnIndexes[j] == -1) {

                    String seColumn = seColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);
                    throw new IllegalArgumentException("Standard error column '" + seColumn + "' not found for '" + variable + "'.");

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

                        double[] result = new double[2 * variableNames.length];

                        for (int j = 0; j < variableNames.length; j++) {

                            String betaString = lineSplit[betaColumnIndexes[j]];

                            double beta = Double.parseDouble(betaString);

                            result[j] = beta;

                            String seString = lineSplit[seColumnIndexes[j]];

                            double se = Double.parseDouble(seString);

                            result[j + variableNames.length] = se;

                        }

                        HashMap<String, HashMap<String, HashMap<String, double[]>>> chromosomeValues = scoringData.get(chromosome);

                        if (chromosomeValues == null) {

                            chromosomeValues = new HashMap<>();
                            scoringData.put(chromosome, chromosomeValues);

                        }

                        HashMap<String, HashMap<String, double[]>> leadVariantValues = chromosomeValues.get(leadVariantId);

                        if (leadVariantValues == null) {

                            leadVariantValues = new HashMap<>(1);
                            chromosomeValues.put(leadVariantId, leadVariantValues);

                        }

                        HashMap<String, double[]> variantValues = leadVariantValues.get(variantId);

                        if (variantValues == null) {

                            variantValues = new HashMap<>(1);
                            leadVariantValues.put(variantId, variantValues);

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
