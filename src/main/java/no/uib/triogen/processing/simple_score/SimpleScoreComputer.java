package no.uib.triogen.processing.simple_score;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
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
     * The type of genotype file.
     */
    private final GenotypesFileType genotypesFileType;
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

    public SimpleScoreComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            ChildToParentMap childToParentMap,
            VariantWeightList variantWeightList,
            File phenotypesFile,
            String[] phenoNames,
            File destinationFile,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.childToParentMap = childToParentMap;
        this.variantWeightList = variantWeightList;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.destinationFile = destinationFile;
        this.logger = logger;

    }

    public void computeScore() {

        logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        HashSet<String> allPhenoColumns = Arrays.stream(phenoNames)
                .collect(
                        Collectors.toCollection(HashSet::new)
                );

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                allPhenoColumns
        );

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Done (" + phenoNames.length + " phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)");

        logger.logMessage("Computing score (geno: " + genotypesFile.getAbsolutePath() + ")");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = GenotypesFileType.getVariantIterator(
                genotypesFile,
                genotypesFileType,
                variantWeightList.variantList,
                0,
                logger
        );

        double[][] scores = new double[childToParentMap.children.length][4];

        GenotypesProvider tempGenotypesProvider;
        while ((tempGenotypesProvider = iterator.next()) != null) {

            GenotypesProvider genotypesProvider = tempGenotypesProvider;
            genotypesProvider.parse(childToParentMap);

            String variantId = genotypesProvider.getVariantID();

            int variantIndex = variantWeightList.variantList.getIndex(variantId);
            String effectAllele = variantWeightList.effectAllele[variantIndex];
            double weight = variantWeightList.weights[variantIndex];

            boolean swapAlleles = effectAllele.equals(genotypesProvider.getRef());

            if (swapAlleles || effectAllele.equals(genotypesProvider.getAlt())) {

                for (int childIndex = 0; childIndex < childToParentMap.children.length; childIndex++) {

                    String childId = childToParentMap.children[childIndex];
                    String motherId = childToParentMap.getMother(childId);
                    String fatherId = childToParentMap.getFather(childId);

                    short[] hs = genotypesProvider.getNAltH(childId, motherId, fatherId);

                    if (!swapAlleles) {

                        for (int haplotypeI = 0; haplotypeI < 4; haplotypeI++) {

                            scores[childIndex][haplotypeI] += weight * hs[haplotypeI];

                        }

                    } else {

                        for (int haplotypeI = 0; haplotypeI < 4; haplotypeI++) {

                            scores[childIndex][haplotypeI] += weight * (1 - hs[haplotypeI]);

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
