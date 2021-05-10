package no.uib.triogen.cmd.prs;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsComputer;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum PrsOptions {

    geno("g", "geno", "The genotypes file.", true, true),
    chromosome("c", "chromosome", "The chromosome name.", true, true),
    phenoFile("p", "phenoFile", "The phenotypes file.", true, true),
    trainingFile("t", "trainingFile", "File listing the summary statistics to use for training.", true, true),
    snpId("s", "snpId", "Name of the variant identifier column in the training file. Default: 'variantId'.", false, true),
    ea("ea", "effectAllele", "Name of the effect allele column in the training file. Default: 'tested_allele'.", false, true),
    model("m", "model", "Names of the model to use. Default: cmf. Available: " + Model.getCommandLineOptions() + ".", false, true),
    variables("v", "variables", "Names of the variables to use, need to be in the same order as specified in the model. Default: c,m,f.", false, true),
    betaPattern("b", "betaPattern", "Pattern for the effect size column. Wildcards: " + PrsComputer.MODEL_WILDCARD + " for model name and " + PrsComputer.VARIABLE_WILDCARD + " for variable name. Default: '" + PrsComputer.MODEL_WILDCARD + ".B" + PrsComputer.VARIABLE_WILDCARD + "'.", false, true),
    childId("id", "childId", "The name of the column containing the child id in the phenotypes file. Default: child_SentrixID.", false, true),
    phenoName("pn", "phenoName", "The names of the phenotype to analyze.", true, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    out("o", "out", "The file where to write the results.", true, true),
    ldThreshold("ld", "ldThreshold", "LD R2 value after which two hits cannot be considered independent. Default: '0.05'.", false, true),
    afThreshold("af", "afThreshold", "Lowest allele frequency considered. 0.001 means that the variants with maf higher or lower than 0.1 % are considered. Default: '0.001'.", false, true);

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
     * @param hasArg has the option an argument
     */
    private PrsOptions(
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

        for (PrsOptions option : values()) {

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
