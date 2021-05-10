package no.uib.triogen.processing.prs;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.covariates.CovariatesHandler;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.Model;

/**
 * This class computes a PRS on trios in a similar way as PRSice for single
 * individuals.
 *
 * @author Marc Vaudel
 */
public class PrsComputer {

    /**
     * Wildcard for the variable name in the summary stats columns.
     */
    public static final String VARIABLE_WILDCARD = "{variable}";
    /**
     * Wildcard for the model name in the summary stats columns.
     */
    public static final String MODEL_WILDCARD = "{model}";
    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The name of the chromosome.
     */
    private final String chromosome;
    /**
     * The file containing the covariates and phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The file to export the result to.
     */
    private final File destinationFile;
    /**
     * The file containing the training summary stats.
     */
    private final File trainingFile;
    /**
     * The column containing the snp id.
     */
    private final String snpIdColumn;
    /**
     * The column containing the effect allele.
     */
    private final String eaColumn;
    /**
     * The trio model to use.
     */
    private final Model model;
    /**
     * The ordered names of the variables used in the model.
     */
    private final String[] variableNames;
    /**
     * The pattern to use to find the effect size column for each variable in
     * the model.
     */
    private final String effectSizeColumnPattern;
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
     * The LD threshold used for pruning.
     */
    private final double ldThreshold;
    /**
     * The minimal allele frequency to consider for thresholding.
     */
    private final double minAf;
    /**
     * The name of the phenotype.
     */
    private final String phenoName;
    /**
     * The name of the covariates to use.
     */
    public final String[] covariates;
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
     * @param chromosome The name of the chromosome.
     * @param phenotypesFile The file containing the covariates and phenotypes.
     * @param inheritanceMap The inheritance map for the given file.
     * @param snpIdColumn The column containing SNP ids.
     * @param trainingFile The file containing the training data.
     * @param model The trio model to use.
     * @param variableNames The names of the variables to include.
     * @param effectSizePattern The pattern to use to find the effect size
     * column for each variable in the model.
     * @param eaColumn The name of the effect allele column.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     * @param childToParentMap The map of trios.
     * @param ldThreshold The LD threshold used for pruning.
     * @param minAf The minimal allele frequency to consider for thresholding.
     * @param phenoName The name of the phenotype.
     * @param covariates The name of the covariates to use.
     * @param destinationFile The file to export the result to.
     * @param logger The logger.
     */
    public PrsComputer(
            File genotypesFile,
            String chromosome,
            File phenotypesFile,
            File destinationFile,
            String snpIdColumn,
            String eaColumn,
            File trainingFile,
            Model model,
            String[] variableNames,
            String effectSizePattern,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy,
            double ldThreshold,
            double minAf,
            String phenoName,
            String[] covariates,
            ChildToParentMap childToParentMap,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.chromosome = chromosome;
        this.phenotypesFile = phenotypesFile;
        this.destinationFile = destinationFile;
        this.trainingFile = trainingFile;
        this.snpIdColumn = snpIdColumn;
        this.eaColumn = eaColumn;
        this.model = model;
        this.variableNames = variableNames;
        this.effectSizeColumnPattern = effectSizePattern;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;
        this.ldThreshold = ldThreshold;
        this.minAf = minAf;
        this.phenoName = phenoName;
        this.covariates = covariates;
        this.childToParentMap = childToParentMap;
        this.logger = logger;

    }

    public void run() throws IOException {

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

        logger.logMessage("Parsing training data from " + trainingFile.getAbsolutePath());

        start = Instant.now().getEpochSecond();

        HashMap<String, Double>[] trainingData = parseTrainingData();

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Parsing training data from " + trainingFile + " done (" + duration + " seconds)");

        logger.logMessage("Parsing covariates and phenotypes from " + phenotypesFile);

        start = Instant.now().getEpochSecond();

        HashSet<String> allPhenoColumns = Arrays.stream(covariates)
                .collect(
                        Collectors.toCollection(HashSet::new)
                );
        allPhenoColumns.add(phenoName);

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                allPhenoColumns
        );

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)");

        logger.logMessage("Adjusting for covariates");

        start = Instant.now().getEpochSecond();

        String[] phenoNames = new String[]{phenoName};
        HashMap<String, TreeSet<String>> covariatesSpecific = new HashMap<>(1);
        covariatesSpecific.put(phenoName, new TreeSet<>());

        CovariatesHandler covariatesHandler = new CovariatesHandler(
                phenoNames,
                phenotypesHandler,
                covariates,
                covariatesSpecific
        );

        phenotypesHandler.adjustForCovariates(phenoNames, covariatesHandler);

        covariatesHandler.trim();

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("    " + phenoName + ": " + covariatesHandler.originalIndexMap.get(phenoName).length + " phenotypes, " + covariatesHandler.covariatesMap.get(phenoName).length + " covariates (including intercept) representing " + covariatesHandler.rankMap.get(phenoName) + " dimensions.");

        logger.logMessage("Done (Adjusted for covariates in " + duration + " seconds)");

        phenotypesHandler.sanityCheck();

    }

    /**
     * For each model variable, parses the training data in a map, variant to
     * effect and p.
     *
     * @return The training data in a map, variant to effect and p, per
     * variable.
     */
    private HashMap<String, Double>[] parseTrainingData() {

        HashMap<String, Double>[] trainingData = new HashMap[variableNames.length];

        for (int j = 0; j < variableNames.length; j++) {

            trainingData[j] = new HashMap<>();

        }

        int snpIdIndex = -1;
        int[] betaColumnIndexes = new int[variableNames.length];

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(trainingFile)) {

            String line = reader.readLine();

            String[] lineSplit = line.split(SEPARATOR);

            for (int i = 0; i < lineSplit.length; i++) {

                if (lineSplit[i].equals(snpIdColumn)) {

                    snpIdIndex = i;

                }

                for (int j = 0; j < variableNames.length; j++) {

                    String variable = variableNames[j];

                    String betaColumn = effectSizeColumnPattern
                            .replace(MODEL_WILDCARD, model.name())
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(betaColumn)) {

                        betaColumnIndexes[j] = i;

                    }
                }
            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(SEPARATOR);

                String snpId = lineSplit[snpIdIndex];

                for (int j = 0; j < variableNames.length; j++) {

                    String betaString = lineSplit[betaColumnIndexes[j]];

                    double beta = Double.parseDouble(betaString);

                    trainingData[j].put(snpId, beta);

                }
            }
        }

        return trainingData;

    }

}
