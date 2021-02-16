package no.uib.triogen.processing.simple_score;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.VariantIterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.simple_score.VariantWeightList;

/**
 * This class computes a simple risk score on trios based on variant weights.
 *
 * @author Marc Vaudel
 */
public class SimpleScoreComputer {

    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The allele inheritance map.
     */
    private final HashMap<Integer, char[]> inheritanceMap;
    /**
     * The default ploidy for mothers.
     */
    private final int defaultMotherPloidy;
    /**
     * The default ploidy for fathers.
     */
    private final int defaultFatherPloidy;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The weight list to use to compute the score.
     */
    private final VariantWeightList variantWeightList;
    /**
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The names of the phenotypes to export.
     */
    private final String[] phenoNames;
    /**
     * The file where to write the results.
     */
    private final File destinationFile;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * The decompressor to use.
     */
    private final ZstdDecompressor decompressor = new ZstdDecompressor();

    public SimpleScoreComputer(
            File genotypesFile,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy,
            ChildToParentMap childToParentMap,
            VariantWeightList variantWeightList,
            File phenotypesFile,
            String[] phenoNames,
            File destinationFile,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;
        this.childToParentMap = childToParentMap;
        this.variantWeightList = variantWeightList;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.destinationFile = destinationFile;
        this.logger = logger;

    }

    public void computeScore() throws IOException {

        logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

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

        logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath());

        start = Instant.now().getEpochSecond();

        HashSet<String> allPhenoColumns = Arrays.stream(phenoNames)
                .collect(
                        Collectors.toCollection(HashSet::new)
                );

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                allPhenoColumns
        );

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (" + phenoNames.length + " phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)");

        logger.logMessage("Computing score (geno: " + genotypesFile.getAbsolutePath() + ")");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = new VariantIterator(
                bgenIndex,
                logger,
                "Simple score in " + genotypesFile.getAbsolutePath()
        );

        double[][] scores = new double[childToParentMap.children.length][4];

        Integer tempIndex;
        while ((tempIndex = iterator.next()) != null) {

            int variantIndex = tempIndex;
            VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

            if (variantInformation.alleles.length > 1) {

                BgenVariantData variantData = bgenFileReader.getVariantData(variantIndex);
                variantData.parse(
                        childToParentMap,
                decompressor
                );

                String variantId = variantInformation.id;
                boolean found = false;

                if (variantWeightList.variantList.contains(variantId)) {

                    found = true;

                }

                if (!found) {

                    variantId = variantInformation.rsId;

                    if (variantWeightList.variantList.contains(variantId)) {

                        found = true;

                    }
                }

                if (found) {

                    int variantWeightIndex = variantWeightList.variantList.getIndex(variantId);
                    String effectAllele = variantWeightList.effectAllele[variantWeightIndex];
                    double weight = variantWeightList.weights[variantWeightIndex];

                    for (int childIndex = 0; childIndex < childToParentMap.children.length; childIndex++) {

                        String childId = childToParentMap.children[childIndex];
                        String motherId = childToParentMap.getMother(childId);
                        String fatherId = childToParentMap.getFather(childId);

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            if (effectAllele.equals(variantInformation.alleles[alleleI])) {

                                double[] haplotypes = variantData.getHaplotypes(childId, motherId, fatherId, alleleI);

                                for (int haplotypeI = 0; haplotypeI < 4; haplotypeI++) {

                                    if (effectAllele.equals(variantInformation.alleles[alleleI])) {

                                        scores[childIndex][haplotypeI] += weight * haplotypes[haplotypeI];

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (Score for " + genotypesFile.getName() + " computed in " + duration + " seconds)");

        logger.logMessage("Exporting score (geno: " + genotypesFile.getAbsolutePath() + ")");

        start = Instant.now().getEpochSecond();

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            StringBuilder line = new StringBuilder();
            line.append("ChildId")
                    .append(IoUtils.SEPARATOR)
                    .append("MaternalTransmitted")
                    .append(IoUtils.SEPARATOR)
                    .append("MaternalNonTransmitted")
                    .append(IoUtils.SEPARATOR)
                    .append("PaternalTransmitted")
                    .append(IoUtils.SEPARATOR)
                    .append("PaternalNonTransmitted");

            for (String phenoName : phenoNames) {

                line.append(IoUtils.SEPARATOR)
                        .append(phenoName);

            }

            writer.writeLine(line.toString());

            for (int childIndex = 0; childIndex < childToParentMap.children.length; childIndex++) {

                line = new StringBuilder();

                String childId = childToParentMap.children[childIndex];

                line.append(childId);

                for (int haplotypeI = 0; haplotypeI < 4; haplotypeI++) {

                    line.append(IoUtils.SEPARATOR)
                            .append(scores[childIndex][haplotypeI]);

                }

                for (String phenoName : phenoNames) {

                    double phenoValue = phenotypesHandler.phenoMap.get(phenoName)[childIndex];

                    line.append(IoUtils.SEPARATOR)
                            .append(phenoValue);

                }

                writer.writeLine(line.toString());

            }
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (Score for " + genotypesFile.getName() + " exported in " + duration + " seconds)");

    }
}
