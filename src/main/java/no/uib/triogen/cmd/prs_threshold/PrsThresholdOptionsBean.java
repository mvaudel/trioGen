package no.uib.triogen.cmd.prs_threshold;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.cmd.prs_score.PrsScoreOptions;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsPruner;
import no.uib.triogen.processing.prs.PrsUtils.ScoringMode;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class PrsThresholdOptionsBean {

    /**
     * The genotypes file path.
     */
    public final String genotypesFile;
    /**
     * The pheno file.
     */
    public final File phenoFile;
    /**
     * The pheno name.
     */
    public final String phenoName;
    /**
     * The trio file.
     */
    public final File trioFile;
    /**
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * The score file.
     */
    public final File scoreFile;
    /**
     * Name of the effect allele column in the training file.
     */
    public Model model = Model.cmf;
    /**
     * Name of the effect allele column in the training file.
     */
    public int nBins = 10;
    /**
     * Names of the variables.
     */
    public String[] variables = new String[]{"c", "m", "f"};
    /**
     * Pattern for the effect size column.
     */
    public String betaPattern = Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD;
    /**
     * Pattern for the standard error column.
     */
    public String sePattern = Model.cmf.name() + ".B.se" + PrsPruner.VARIABLE_WILDCARD;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * The highest p-value to consider.
     */
    public double pValueThreshold = 1e-6;
    /**
     * The allele frequency threshold.
     */
    public double afThreshold = Double.NaN;
    /**
     * Column containing the allele frequency.
     */
    public String afColumn = null;
    /**
     * The quantile to use for beta estimation.
     */
    public double betaQuantile = 0.025;
    /**
     * The scoring mode
     */
    public ScoringMode scoringMode = ScoringMode.weighted;
    /**
     * The folder where to store the bgen index files.
     */
    public File bgenIndexFolder;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public PrsThresholdOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsThresholdOptions option : PrsThresholdOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        genotypesFile = aLine.getOptionValue(PrsThresholdOptions.geno.opt);

        File genotypesFolder = (new File(genotypesFile)).getParentFile();

        if (!genotypesFolder.exists()) {

            throw new IllegalArgumentException("Folder containing the genotypes file (" + genotypesFolder + ") not found.");

        } else {

            bgenIndexFolder = genotypesFolder;

        }

        // The phenotypes file
        String filePath = aLine.getOptionValue(PrsThresholdOptions.phenoFile.opt);
        
        phenoFile = new File(filePath);

        if (!phenoFile.exists()) {

            throw new IllegalArgumentException("Phenotypes file (" + genotypesFolder + ") not found.");

        }

        // The phenotype name
        phenoName = aLine.getOptionValue(PrsThresholdOptions.phenoName.opt);
        

        // The genotypes file
        if (aLine.hasOption(PrsThresholdOptions.bgenIndexFolder.opt)) {

            filePath = aLine.getOptionValue(PrsThresholdOptions.bgenIndexFolder.opt);

            bgenIndexFolder = new File(filePath);

            if (!bgenIndexFolder.exists()) {

                throw new IllegalArgumentException("Bgen index folder (" + bgenIndexFolder + ") not found.");

            }
        }

        // the trio file
        filePath = aLine.getOptionValue(PrsThresholdOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The variant ids
        if (aLine.hasOption(PrsThresholdOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(PrsThresholdOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The summary stats file
        filePath = aLine.getOptionValue(PrsThresholdOptions.scoreFile.opt);

        scoreFile = new File(filePath);

        if (!scoreFile.exists()) {

            throw new IllegalArgumentException("Training file (" + scoreFile + ") not found.");

        }

        // The scoring mode
        if (aLine.hasOption(PrsThresholdOptions.model.opt)) {

            String option = aLine.getOptionValue(PrsThresholdOptions.model.opt);

            model = Model.valueOf(option);

        }

        // The variables
        if (aLine.hasOption(PrsThresholdOptions.variables.opt)) {

            String option = aLine.getOptionValue(PrsThresholdOptions.variables.opt);

            variables = option.split(",");

            if (variables.length != model.betaNames.length) {

                String modelVariables = Arrays.stream(variables).collect(Collectors.joining(","));

                throw new IllegalArgumentException("Found " + variables.length + " variables (" + option + ") where " + model.betaNames.length + " (" + modelVariables + ") expected.");

            }
        }

        // The beta pattern
        if (aLine.hasOption(PrsThresholdOptions.betaPattern.opt)) {

            betaPattern = aLine.getOptionValue(PrsThresholdOptions.betaPattern.opt);

        }

        // The se pattern
        if (aLine.hasOption(PrsThresholdOptions.sePattern.opt)) {

            sePattern = aLine.getOptionValue(PrsThresholdOptions.sePattern.opt);

        }

        // The af threshold
        if (aLine.hasOption(PrsScoreOptions.afThreshold.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.afThreshold.opt);

            try {

                afThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("The value for allele frequency threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(afThreshold) || Double.isInfinite(afThreshold)) {

                throw new IllegalArgumentException("The allele frequency threshold (" + option + ") could not be parsed as a number.");

            }
            if (afThreshold <= 0 || afThreshold >= 0.5) {

                throw new IllegalArgumentException("The allele frequency threshold (" + option + ") should be higher than 0 and lower than 0.5.");

            }
            
            if (!aLine.hasOption(PrsScoreOptions.afColumn.opt)) {

                throw new IllegalArgumentException("No column provided for allele frequency.");
                
            }
            
            afColumn = aLine.getOptionValue(PrsScoreOptions.afColumn.opt);
            
        }

        // The output file
        filePath = aLine.getOptionValue(PrsThresholdOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
