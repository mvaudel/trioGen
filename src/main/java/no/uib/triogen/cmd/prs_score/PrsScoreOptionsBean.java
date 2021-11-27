package no.uib.triogen.cmd.prs_score;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.cmd.prs_threshold.PrsThresholdOptions;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsPruner;
import no.uib.triogen.processing.prs.PrsUtils.ScoringMode;
import no.uib.triogen.utils.cli.CliUtils;
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
    public PrsScoreOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsScoreOptions option : PrsScoreOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        genotypesFile = CliUtils.getOptionValue(aLine, PrsScoreOptions.geno);

        File genotypesFolder = (new File(genotypesFile)).getParentFile();

        if (!genotypesFolder.exists()) {

            throw new IllegalArgumentException("Folder containing the genotypes file (" + genotypesFolder + ") not found.");

        } else {

            bgenIndexFolder = genotypesFolder;

        }

        // The genotypes file
        if (CliUtils.hasOption(aLine, PrsScoreOptions.bgenIndexFolder)) {

            String filePath = CliUtils.getOptionValue(aLine, PrsScoreOptions.bgenIndexFolder);

            bgenIndexFolder = new File(filePath);

            if (!bgenIndexFolder.exists()) {

                throw new IllegalArgumentException("Bgen index folder (" + bgenIndexFolder + ") not found.");

            }
        }

        // the trio file
        String filePath = CliUtils.getOptionValue(aLine, PrsScoreOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The variant ids
        if (CliUtils.hasOption(aLine, PrsScoreOptions.variantId)) {

            filePath = CliUtils.getOptionValue(aLine, PrsScoreOptions.variantId);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The summary stats file
        filePath = CliUtils.getOptionValue(aLine, PrsScoreOptions.scoreFile);

        scoreFile = new File(filePath);

        if (!scoreFile.exists()) {

            throw new IllegalArgumentException("Training file (" + scoreFile + ") not found.");

        }

        // The scoring mode
        if (CliUtils.hasOption(aLine, PrsScoreOptions.scoringMode)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.model);

            model = Model.valueOf(option);

        }

        // The variables
        if (CliUtils.hasOption(aLine, PrsScoreOptions.variables)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.variables);

            variables = option.split(",");

            if (variables.length != model.betaNames.length) {

                String modelVariables = Arrays.stream(variables).collect(Collectors.joining(","));

                throw new IllegalArgumentException("Found " + variables.length + " variables (" + option + ") where " + model.betaNames.length + " (" + modelVariables + ") expected.");

            }
        }

        // The model
        if (CliUtils.hasOption(aLine, PrsScoreOptions.scoringMode)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.scoringMode);

            int selectedOption;
            try {

                selectedOption = Integer.parseInt(option);

            } catch (Exception e) {

                throw new IllegalArgumentException("The value provided for scoring mode ('" + option + "') could not be parsed as a number.");

            }

            boolean found = false;

            for (ScoringMode scoringModeOption : ScoringMode.values()) {

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
        if (CliUtils.hasOption(aLine, PrsScoreOptions.betaPattern)) {

            betaPattern = CliUtils.getOptionValue(aLine, PrsScoreOptions.betaPattern);

        }

        // The se pattern
        if (CliUtils.hasOption(aLine, PrsThresholdOptions.sePattern)) {

            sePattern = CliUtils.getOptionValue(aLine, PrsThresholdOptions.sePattern);

        }

        // The p-value threshold
        if (CliUtils.hasOption(aLine, PrsScoreOptions.pValueThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.pValueThreshold);

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

        // The af threshold
        if (CliUtils.hasOption(aLine, PrsScoreOptions.afThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.afThreshold);

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
            
            if (!CliUtils.hasOption(aLine, PrsScoreOptions.afColumn)) {

                throw new IllegalArgumentException("No column provided for allele frequency.");
                
            }
            
            afColumn = CliUtils.getOptionValue(aLine, PrsScoreOptions.afColumn);
            
        }

        // The beta quantile
        if (CliUtils.hasOption(aLine, PrsScoreOptions.betaQuantile)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.betaQuantile);

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
        filePath = CliUtils.getOptionValue(aLine, PrsScoreOptions.out);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
