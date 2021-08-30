package no.uib.triogen.processing.prs;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.utils.SimpleSemaphore;
import static no.uib.triogen.utils.Utils.CHROMOSOME_WILDCARD;

/**
 * This class computes a score for the different PRS classes.
 *
 * @author Marc Vaudel
 */
public class PrsComputer {

    /**
     * Returns the scores for the given trios.
     * 
     * @param genotypesFilePath The path to the genotype files.
     * @param scoringData The scoring data.
     * @param childToParentMap The child to parent map.
     * @param model The model to use.
     * @param variableNames The name of the variables.
     * @param seScaling The se scaling factor.
     * @param bgenIndexFolder The folder where to store the bgen index files.
     * @param semaphore The semaphore to use for the edition of the score map.
     * @param logger The logger to use to display feedback and errors.
     * 
     * @return The scores for this model variables indexed by child id.
     */
    public static HashMap<String, double[]> getScores(
            String genotypesFilePath,
            HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData,
            ChildToParentMap childToParentMap,
            Model model,
            String[] variableNames,
            double seScaling,
            File bgenIndexFolder,
            SimpleSemaphore semaphore,
            SimpleCliLogger logger
    ) {

        HashMap<String, double[]> scores = new HashMap<>(childToParentMap.children.length);

        for (String childId : childToParentMap.children) {

            scores.put(childId, new double[variableNames.length]);

        }

        scoringData.entrySet()
                .stream()
                .parallel()
                .forEach(
                        entry -> processChromosome(
                                childToParentMap,
                                entry.getKey(),
                                genotypesFilePath,
                                model,
                                variableNames,
                                seScaling,
                                entry.getValue(),
                                scores,
                                bgenIndexFolder,
                                semaphore,
                                logger
                        )
                );

        return scores;

    }

    /**
     * Processes the given chromosome.
     * 
     * @param genotypesFilePath The path to the genotype files.
     * @param chromosome The chromosome to process.
     * @param chromosomeScoringData The scoring data for this chromosome.
     * @param scores The map where to store the scores.
     * @param childToParentMap The child to parent map.
     * @param model The model to use.
     * @param variableNames The name of the variables.
     * @param seScaling The se scaling factor.
     * @param bgenIndexFolder The folder where to store the bgen index files.
     * @param semaphore The semaphore to use for the edition of the score map.
     * @param logger The logger to use to display feedback and errors.
     */
    private static void processChromosome(
            ChildToParentMap childToParentMap,
            String chromosome,
            String genotypesFilePath,
            Model model,
            String[] variableNames,
            double seScaling,
            HashMap<String, HashMap<String, HashMap<String, double[]>>> chromosomeScoringData,
            HashMap<String, double[]> scores,
            File bgenIndexFolder,
            SimpleSemaphore semaphore,
            SimpleCliLogger logger
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

            BgenIndex bgenIndex = BgenIndex.getBgenIndex(
                    genotypesFile,
                    BgenIndex.getDefaultIndexFile(
                            bgenIndexFolder,
                            genotypesFile
                    )
            );
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

            HashSet<String> variantsScore = new HashSet<>();
            int nLociScore = chromosomeScoringData.size();

            start = Instant.now().getEpochSecond();

            HashMap<String, ArrayList<String>> scoringToLeadVariantMap = new HashMap<>();

            for (Map.Entry<String, HashMap<String, HashMap<String, double[]>>> entry : chromosomeScoringData.entrySet()) {

                String leadVariant = entry.getKey();

                for (String scoringVariant : entry.getValue().keySet()) {

                    variantsScore.add(scoringVariant);

                    ArrayList<String> leadVariants = scoringToLeadVariantMap.get(scoringVariant);

                    if (leadVariants == null) {

                        leadVariants = new ArrayList<>(1);
                        scoringToLeadVariantMap.put(scoringVariant, leadVariants);

                    }

                    leadVariants.add(leadVariant);

                }
            }

            int nVariantsScore = variantsScore.size();

            HashMap<String, HashMap<String, Integer>> leadToScoringVariantIndexMap = new HashMap<>();
            int nVariants = 0;

            for (int bgenI = 0; bgenI < bgenIndex.variantIdArray.length; bgenI++) {

                String bgenVariantId = bgenIndex.variantIdArray[bgenI];

                ArrayList<String> leadVariants = scoringToLeadVariantMap.get(bgenVariantId);

                if (leadVariants != null) {

                    for (String leadVariant : leadVariants) {

                        HashMap<String, Integer> scoringVariants = leadToScoringVariantIndexMap.get(leadVariant);

                        if (scoringVariants == null) {

                            scoringVariants = new HashMap<>(1);
                            leadToScoringVariantIndexMap.put(leadVariant, scoringVariants);

                        }

                        scoringVariants.put(bgenVariantId, bgenI);

                    }

                    nVariants++;

                }
            }

            int nLociFound = leadToScoringVariantIndexMap.size();

            end = Instant.now().getEpochSecond();
            duration = end - start;

            logger.logMessage("Gathering genotyped variants from " + chromosome + " done (" + duration + " seconds), " + nVariants + " variants of " + nVariantsScore + ", " + nLociFound + " loci of " + nLociScore + ".");

            logger.logMessage("Scoring chromosome " + chromosome);

            start = Instant.now().getEpochSecond();

            leadToScoringVariantIndexMap.entrySet()
                    .stream()
                    .parallel()
                    .forEach(
                            entry -> {

                                String leadVariant = entry.getKey();
                                HashMap<String, Integer> indexMap = entry.getValue();
                                double weight = 1.0 / ((double) indexMap.size());

                                for (Map.Entry<String, Integer> variantEntry : indexMap.entrySet()) {

                                    String scoringVariant = variantEntry.getKey();
                                    int bgenVariantIndex = variantEntry.getValue();

                                    HashMap<String, double[]> alleles = chromosomeScoringData.get(leadVariant).get(scoringVariant);

                                    for (Map.Entry<String, double[]> allelesEntry : alleles.entrySet()) {

                                        String allele = allelesEntry.getKey();
                                        double[] scoringDetails = allelesEntry.getValue();

                                        double[] betaContributions = new double[variableNames.length];

                                        boolean nonNull = false;

                                        for (int variableI = 0; variableI < variableNames.length; variableI++) {

                                            double beta = scoringDetails[variableI];
                                            double se = scoringDetails[variableI + variableNames.length];

                                            double betaEstimate = 0.0;

                                            if (beta > 0.0) {

                                                betaEstimate = beta + seScaling * se;

                                                if (betaEstimate <= 0.0) {

                                                    betaEstimate = 0.0;

                                                } else {

                                                    nonNull = true;

                                                }
                                            } else if (beta < 0.0) {

                                                betaEstimate = beta - seScaling * se;

                                                if (betaEstimate >= 0.0) {

                                                    betaEstimate = 0.0;

                                                } else {

                                                    nonNull = true;

                                                }
                                            }

                                            betaContributions[variableI] = betaEstimate;

                                        }

                                        if (nonNull) {

                                            score(
                                                    childToParentMap,
                                                    model,
                                                    bgenVariantIndex,
                                                    allele,
                                                    weight,
                                                    betaContributions,
                                                    bgenFileReader,
                                                    bgenIndex,
                                                    semaphore,
                                                    scores
                                            );
                                        }
                                    }
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


    /**
     * Computes the score contribution for the given variant.
     * 
     * @param scores The map where to store the scores.
     * @param childToParentMap The child to parent map.
     * @param model The model to use.
     * @param bgenVariantIndex The index of the variant in the bgen file.
     * @param effectAllele The effect allele.
     * @param weight The weight for this variant.
     * @param betaContributions The beta contributions for the different variables.
     * @param bgenFileReader The bgen file reader.
     * @param bgenIndex The index of the bgen file.
     * @param semaphore The semaphore to use for the edition of the score map.
     */
    private static void score(
            ChildToParentMap childToParentMap,
            Model model,
            int bgenVariantIndex,
            String effectAllele,
            double weight,
            double[] betaContributions,
            BgenFileReader bgenFileReader,
            BgenIndex bgenIndex,
            SimpleSemaphore semaphore,
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

                double betaEstimate = betaContributions[variableI];

                double variantContribution = xValue * betaEstimate * weight;

                childScores[variableI] = childScores[variableI] + variantContribution;

            }

            semaphore.release();

        }
    }
}
