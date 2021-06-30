package no.uib.triogen.cmd.prs_score;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsScorer;
import no.uib.triogen.processing.prs.PrsTrainer;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum PrsScoreOptions {

    geno("g", "geno", "The genotypes file. Wildcard: " + PrsScorer.CHROMOSOME_WILDCARD + " for variable name.", true, true),
    variantId("vi", "variantId", "File listing the variants to include in the analysis. Default: process all variants in the genotypes file.", false, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    scoreFile("s", "scoreFile", "Score details as obtained from the PrsTrain command.", true, true),
    model("m", "model", "Names of the model to use. Default: " + Model.cmf.name() + ". Available: " + Model.getCommandLineOptions() + ".", false, true),
    variables("v", "variables", "Names of the variables to use, need to be in the same order as specified in the model. Default: c,m,f.", false, true),
    betaPattern("b", "betaPattern", "Pattern for the effect size column. Wildcard: " + PrsTrainer.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsTrainer.VARIABLE_WILDCARD + "'.", false, true),
    sePattern("s", "sePattern", "Pattern for the standard error column. Wildcard: " + PrsTrainer.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsTrainer.VARIABLE_WILDCARD + ".se'.", false, true),
    scoringMode("sm", "scoringMode", "The scoring mode. 0: top hit per locus; 1: weighted average. Default: 1.", false, true),
    out("o", "out", "The file where to write the results.", true, true),
    pValueThreshold("pv", "pValueThreshold", "The highest p-value to consider. Default: '1e-6'.", false, true),
    betaQuantile("bq", "betaQuantile", "The probability quantile for the beta estimation. Default: '0.025'.", false, true);

    /**
     * The short option.
     */
    public final String opt;
    /**
     * The long option.
     */
    public final String longOpt;
    /**
     * Explanation for the CLI option.
     */
    public final String description;
    /**
     * Boolean indicating whether the option is mandatory.
     */
    public final boolean mandatory;
    /**
     * Boolean indicating whether the option has an argument.
     */
    public final boolean hasArg;

    /**
     * Private constructor managing the various variables for the enum
     * instances.
     *
     * @param opt the sort option
     * @param longOpt the long option
     * @param description the description
     * @param mandatory is the option mandatory
     * 
     * 
     * @param hasArg has the option an argument
     */
    private PrsScoreOptions(
            String opt, 
            String longOpt, 
            String description, 
            boolean mandatory, 
            boolean hasArg
    ) {
        this.opt = opt;
        this.longOpt = longOpt;
        this.description = description;
        this.mandatory = mandatory;
        this.hasArg = hasArg;
    }

    /**
     * Creates the options for the command line interface based on the possible
     * values.
     *
     * @param options the apache options object
     */
    public static void createOptionsCLI(
            Options options
    ) {

        for (PrsScoreOptions option : values()) {

            options.addOption(option.opt, option.longOpt, option.hasArg, option.description);

        }
    }

    /**
     * Returns the options as a string.
     *
     * @return the options as a string
     */
    public static String getOptionsAsString() {

        final StringBuilder output = new StringBuilder();
        String formatter = "%-35s";

        output.append("General Options:");
        output.append(LINE_SEPARATOR)
                .append(LINE_SEPARATOR);
        
        output.append("-").append(String.format(formatter, "h (--help)")).append(" ").append("Shows a brief help message.").append(LINE_SEPARATOR);
        output.append("-").append(String.format(formatter, "v (--version)")).append(" ").append("Shows the version of the tool.").append(LINE_SEPARATOR);

        output.append(LINE_SEPARATOR)
                .append(LINE_SEPARATOR);
        output.append("Mandatory Options:");
        output.append(LINE_SEPARATOR)
                .append(LINE_SEPARATOR);

        Arrays.stream(values())
                .filter(option -> option.mandatory)
                .forEach(option -> output.append("-").append(String.format(formatter, option.opt + " (--" + option.longOpt + ")")).append(" ").append(option.description).append(LINE_SEPARATOR));

        output.append(LINE_SEPARATOR)
                .append(LINE_SEPARATOR);
        output.append("Additional Options:");
        output.append(LINE_SEPARATOR)
                .append(LINE_SEPARATOR);

        Arrays.stream(values())
                .filter(option -> !option.mandatory)
                .forEach(option -> output.append("-").append(String.format(formatter, option.opt + " (--" + option.longOpt + ")")).append(" ").append(option.description).append(LINE_SEPARATOR));

        return output.toString();
    }
}
