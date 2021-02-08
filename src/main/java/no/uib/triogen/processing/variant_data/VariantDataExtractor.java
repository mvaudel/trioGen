package no.uib.triogen.processing.variant_data;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.VariantIterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.covariates.CovariatesHandler;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class extracts data for a given sets of variants.
 *
 * @author Marc Vaudel
 */
public class VariantDataExtractor {

    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The inheritance map.
     */
    private final HashMap<Integer, char[]> inheritanceMap;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The names of the phenotypes to use.
     */
    private final String[] phenoNames;
    /**
     * Names of the covariates to use for all the phenotypes.
     */
    public final String[] covariatesGeneral;
    /**
     * Map of the covariates to use for specific phenotypes.
     */
    public final HashMap<String, TreeSet<String>> covariatesSpecific;
    /**
     * The file to export the result to.
     */
    private final File destinationFile;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * Constructor.
     *
     * @param genotypesFile The file containing the genotypes.
     * @param inheritanceMap The inheritance map for the given file.
     * @param variantList The variants to process.
     * @param childToParentMap The map of trios.
     * @param phenotypesFile The file containing the phenotypes.
     * @param phenoNames The names of the phenotypes to use.
     * @param covariatesGeneral The names of the general covariates to use for
     * all phenotypes.
     * @param covariatesSpecific The names of the covariates specific to
     * phenotypes.
     * @param destinationFile The file to export the result to.
     * @param logger The logger.
     */
    public VariantDataExtractor(
            File genotypesFile,
            HashMap<Integer, char[]> inheritanceMap,
            VariantList variantList,
            ChildToParentMap childToParentMap,
            File phenotypesFile,
            String[] phenoNames,
            String[] covariatesGeneral,
            HashMap<String, TreeSet<String>> covariatesSpecific,
            File destinationFile,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.inheritanceMap = inheritanceMap;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.covariatesGeneral = covariatesGeneral;
        this.covariatesSpecific = covariatesSpecific;
        this.destinationFile = destinationFile;
        this.logger = logger;

    }

    /**
     * Goes through the list of provided variants and extracts the values of the
     * adjusted phenotypes.
     * 
     * @throws java.io.IOException Exception thrown if an error occurs while reading or writing a file.
     */
    public void run() throws IOException {

        logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        BgenIndex bgenIndex = BgenIndex.getBgenIndex(genotypesFile);
        BgenFileReader bgenFileReader = new BgenFileReader(genotypesFile, bgenIndex, inheritanceMap);

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing " + genotypesFile + " done (" + duration + " seconds)");

        logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath());

        start = Instant.now().getEpochSecond();

        HashSet<String> allPhenoColumns = Stream.concat(
                Arrays.stream(phenoNames),
                Stream.concat(
                        Arrays.stream(covariatesGeneral),
                        Stream.concat(
                                covariatesSpecific.keySet().stream(),
                                covariatesSpecific.values().stream()
                                        .flatMap(
                                                subCovariates -> subCovariates.stream()
                                        )
                        )
                )
        ).collect(
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

        logger.logMessage("Adjusting for covariates");

        start = Instant.now().getEpochSecond();

        CovariatesHandler covariatesHandler = new CovariatesHandler(
                phenoNames,
                phenotypesHandler,
                covariatesGeneral,
                covariatesSpecific
        );

        phenotypesHandler.adjustForCovariates(phenoNames, covariatesHandler);

        covariatesHandler.trim();

        end = Instant.now().getEpochSecond();
        duration = end - start;

        Arrays.stream(phenoNames)
                .sorted()
                .forEach(
                        phenoName -> logger.logMessage("    " + phenoName + ": " + covariatesHandler.originalIndexMap.get(phenoName).length + " phenotypes, " + covariatesHandler.covariatesMap.get(phenoName).length + " covariates (including intercept) representing " + covariatesHandler.rankMap.get(phenoName) + " dimensions.")
                );

        logger.logMessage("Done (Adjusted for covariates in " + duration + " seconds)");

        phenotypesHandler.sanityCheck();

        logger.logMessage("Exctracting data (geno: " + genotypesFile.getAbsolutePath() + variantList.variantId.length + " variants, pheno: " + phenotypesFile.getAbsolutePath() + " " + phenoNames.length + " phenotypes)");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = new VariantIterator(
                bgenIndex,
                logger,
                "Mendelian error check in " + genotypesFile.getAbsolutePath()
        );

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            writer.writeLine(
                    "contig",
                    "position",
                    "variantId",
                    "rsId",
                    "allele",
                    "otherAllele",
                    "motherNonTransmitted",
                    "motherTransmitted",
                    "fatherTransmitted",
                    "fatherNonTransmitted",
                    String.join(
                            SEPARATOR,
                            phenoNames
                    )
            );

            Integer tempIndex;
            while ((tempIndex = iterator.next()) != null) {

                int variantIndex = tempIndex;
                VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

                if (variantInformation.alleles.length > 1) {

                    BgenVariantData variantData = bgenFileReader.getVariantData(variantIndex);
                    variantData.parse(childToParentMap);

                    for (int i = 0; i < childToParentMap.children.length; i++) {

                        int trioIndex = i;
                        String childId = childToParentMap.children[trioIndex];
                        String motherId = childToParentMap.getMother(childId);
                        String fatherId = childToParentMap.getFather(childId);

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            double[] haplotypes = variantData.getHaplotypes(childId, motherId, fatherId, alleleI);

                            StringBuilder phenoValues = new StringBuilder();

                            for (String phenoName : phenoNames) {

                                if (phenoValues.length() > 0) {

                                    phenoValues.append(SEPARATOR);

                                }

                                int adjustedIndex = covariatesHandler.adjustedIndexMap.get(phenoName)[trioIndex];

                                if (adjustedIndex != -1) {

                                    double[] adjustedValues = phenotypesHandler.phenoMap.get(phenoName);
                                    double phenoValue = adjustedValues[adjustedIndex];

                                    phenoValues.append(phenoValue);

                                } else {

                                    phenoValues.append(Double.NaN);

                                }
                            }

                            writer.writeLine(
                                    variantInformation.contig,
                                    Integer.toString(variantInformation.position),
                                    variantInformation.id,
                                    variantInformation.rsId,
                                    variantInformation.alleles[alleleI],
                                    variantInformation.getOtherAllele(alleleI),
                                    Double.toString(haplotypes[0]),
                                    Double.toString(haplotypes[1]),
                                    Double.toString(haplotypes[2]),
                                    Double.toString(haplotypes[3]),
                                    phenoValues.toString()
                            );
                        }
                    }
                }
            }
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (Mendelian error check in " + duration + " seconds)");
        
    }
}
