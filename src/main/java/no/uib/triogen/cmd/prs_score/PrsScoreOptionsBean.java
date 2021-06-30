package no.uib.triogen.cmd.prs_score;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.cmd.prs_train.PrsTrainOptions;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsScorer;
import no.uib.triogen.processing.prs.PrsTrainer;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class PrsScoreOptionsBean {

    /**
     * The genotypes file path.
     */
    public final String genotypesFile;
    /**
     * the trio file.
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
     * Names of the variables.
     */
    public String[] variables = new String[]{"c", "m", "f"};
    /**
     * Pattern for the effect size column.
     */
    public String betaPattern = Model.cmf.name() + ".B" + PrsTrainer.VARIABLE_WILDCARD;
    /**
     * Pattern for the standard error column.
     */
    public String sePattern = Model.cmf.name() + ".B.se" + PrsTrainer.VARIABLE_WILDCARD;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * The highest p-value to consider.
     */
    public double pValueThreshold = 1e-6;
    /**
     * The quantile to use for beta estimation.
     */
    public double betaQuantile = 0.025;
    /**
     * The scoring mode
     */
    public PrsScorer.ScoringMode scoringMode = PrsScorer.ScoringMode.weighted;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public PrsScoreOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsScoreOptions option : PrsScoreOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        genotypesFile = aLine.getOptionValue(PrsScoreOptions.geno.opt);
        
        File genotypesFolder = (new File(genotypesFile)).getParentFile();

        if (!genotypesFolder.exists()) {

            throw new IllegalArgumentException("Folder containing the genotypes file (" + genotypesFolder + ") not found.");

        }

        // the trio file
        String filePath = aLine.getOptionValue(PrsScoreOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The variant ids
        if (aLine.hasOption(PrsScoreOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(PrsScoreOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The summary stats file
        filePath = aLine.getOptionValue(PrsScoreOptions.scoreFile.opt);

        scoreFile = new File(filePath);

        if (!scoreFile.exists()) {

            throw new IllegalArgumentException("Training file (" + scoreFile + ") not found.");

        }
        
        // The scoring mode
        if (aLine.hasOption(PrsScoreOptions.scoringMode.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.model.opt);

            model = Model.valueOf(option);

        }

        // The variables
        if (aLine.hasOption(PrsScoreOptions.variables.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.variables.opt);

            variables = option.split(",");

            if (variables.length != model.betaNames.length) {

                String modelVariables = Arrays.stream(variables).collect(Collectors.joining(","));

                throw new IllegalArgumentException("Found " + variables.length + " variables (" + option + ") where " + model.betaNames.length + " (" + modelVariables + ") expected.");

            }
        }

        // The model
        if (aLine.hasOption(PrsScoreOptions.scoringMode.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.scoringMode.opt);
            
            int selectedOption;
            try {
            
             selectedOption = Integer.parseInt(option);
            
            } catch (Exception e) {
                
                throw new IllegalArgumentException("The value provided for scoring mode ('" + option + "') could not be parsed as a number.");
                
            }
            
            boolean found = false;
            
            for (PrsScorer.ScoringMode scoringModeOption : PrsScorer.ScoringMode.values()) {
                
                if (scoringModeOption.index == selectedOption) {
                    
                    scoringMode = scoringModeOption;
                    found = true;
                    
                    break;
                    
                }
            }
            
            if (!found) {
                
                throw new IllegalArgumentException("The value provided for scoring mode ('" + option + "') does not correspond to a scoring mode.");
                
            }

        }

        // The beta pattern
        if (aLine.hasOption(PrsScoreOptions.betaPattern.opt)) {

            betaPattern = aLine.getOptionValue(PrsScoreOptions.betaPattern.opt);

        }

        // The se pattern
        if (aLine.hasOption(PrsTrainOptions.sePattern.opt)) {

            sePattern = aLine.getOptionValue(PrsTrainOptions.sePattern.opt);

        }

        // The p-value threshold
        if (aLine.hasOption(PrsScoreOptions.pValueThreshold.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.pValueThreshold.opt);

            try {

                pValueThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("The value for p-value threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(pValueThreshold) || Double.isInfinite(pValueThreshold)) {

                throw new IllegalArgumentException("The p-value threshold (" + option + ") could not be parsed as a number.");

            }
            if (pValueThreshold <= 0 || pValueThreshold > 1) {

                throw new IllegalArgumentException("The p-value threshold (" + option + ") should be higher than 0 and lower than 1.");

            }
        }

        // The beta quantile
        if (aLine.hasOption(PrsScoreOptions.betaQuantile.opt)) {

            String option = aLine.getOptionValue(PrsScoreOptions.betaQuantile.opt);

            try {

                betaQuantile = Double.parseDouble(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("The value for beta quantile (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(betaQuantile) || Double.isInfinite(betaQuantile)) {

                throw new IllegalArgumentException("The value for beta quantile (" + option + ") could not be parsed as a number.");

            }
            if (betaQuantile <= 0 || betaQuantile > 0.5) {

                throw new IllegalArgumentException("The value for beta quantile (" + option + ") should be higher than 0 and lower or equal to 0.5.");

            }
        }

        // The output file
        filePath = aLine.getOptionValue(PrsScoreOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
