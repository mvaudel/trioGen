package no.uib.triogen.processing.prs;

import java.io.File;
import java.util.HashMap;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.Model;

/**
 * This class computes a PRS on trios in a similar way as PRSice for single
 * individuals.
 *
 * @author Marc Vaudel
 */
public class PrsComputer {

    public static final String VARIABLE_WILDCARD = "{variable}";
    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
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
     * The pattern to use to find the p-value column for each variable in the
     * model.
     */
    private final String pValueColumnPattern;
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
     * The maximal allele frequency to consider for thresholding.
     */
    private final double maxAf;
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
     * @param inheritanceMap The inheritance map for the given file.
     * @param model The trio model to use.
     * @param effectSizePattern The pattern to use to find the effect size
     * column for each variable in the model.
     * @param pValuePattern The pattern to use to find the p-value column for
     * each variable in the model.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     * @param childToParentMap The map of trios.
     * @param ldThreshold The LD threshold used for pruning.
     * @param minAf The minimal allele frequency to consider for thresholding.
     * @param maxAf The maximal allele frequency to consider for thresholding.
     * @param covariates The name of the covariates to use.
     * @param destinationFile The file to export the result to.
     * @param logger The logger.
     */
    public PrsComputer(
            File genotypesFile,
            File destinationFile,
            String snpIdColumn,
            File trainingFile,
            Model model,
            String[] variableNames,
            String effectSizePattern,
            String pValuePattern,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy,
            double ldThreshold,
            double minAf,
            double maxAf,
            String[] covariates,
            ChildToParentMap childToParentMap,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.destinationFile = destinationFile;
        this.trainingFile = trainingFile;
        this.snpIdColumn = snpIdColumn;
        this.model = model;
        this.variableNames = variableNames;
        this.effectSizeColumnPattern = effectSizePattern;
        this.pValueColumnPattern = pValuePattern;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;
        this.ldThreshold = ldThreshold;
        this.minAf = minAf;
        this.maxAf = maxAf;
        this.covariates = covariates;
        this.childToParentMap = childToParentMap;
        this.logger = logger;

    }

    public void run() {

    }

    /**
     * For each model variable, parses the training data in a map, variant to
     * effect and p.
     *
     * @return The training data in a map, variant to effect and p, per
     * variable.
     */
    private HashMap<String, double[]>[] parseTrainingData() {

        HashMap<String, double[]>[] trainingData = new HashMap[variableNames.length];

                for (int j = 0; j < variableNames.length; j++) {
                    
                    trainingData[j] = new HashMap<>();
                    
                }

        int snpIdIndex = -1;
        int[] betaColumnIndexes = new int[variableNames.length];
        int[] pColumnIndexes = new int[variableNames.length];

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(genotypesFile)) {

            String line = reader.readLine();

            String[] lineSplit = line.split(SEPARATOR);

            for (int i = 0; i < lineSplit.length; i++) {

                if (lineSplit[i].equals(snpIdColumn)) {

                    snpIdIndex = i;

                }

                for (int j = 0; j < variableNames.length; j++) {

                    String variable = variableNames[j];

                    String betaColumn = effectSizeColumnPattern.replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(betaColumn)) {

                        betaColumnIndexes[j] = i;

                    }

                    String pColumn = pValueColumnPattern.replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(pColumn)) {

                        pColumnIndexes[j] = i;

                    }
                }
            }
            
            while ((line = reader.readLine()) != null) {
                
                lineSplit = line.split(SEPARATOR);
                
                String snpId = lineSplit[snpIdIndex];

                for (int j = 0; j < variableNames.length; j++) {
                    
                    String betaString = lineSplit[betaColumnIndexes[j]];
                    
                    double beta = Double.parseDouble(betaString);
                    
                    String pString = lineSplit[pColumnIndexes[j]];
                    
                    double p = Double.parseDouble(pString);
                    
                    double[] variableValues = new double[]{beta, p};
                    
                    trainingData[j].put(snpId, variableValues);
                    
                }
            }
        }
        
        return trainingData;

    }

}
