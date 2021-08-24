package no.uib.triogen.cmd.ld_pruning;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.utils.Utils;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum LdPruningOptions {

    results("res", "results", "The results file to prune.", true, true),
    ldMatrix("l", "ldMatrix", "The ld matrix file as generated using the LdMatrix command. If LD matrix files are computed per contig, replace the chromosome name with '" + Utils.CHROMOSOME_WILDCARD + "'.", true, true),
    minR2("r", "minR2", "The minimal ld r2 to consider two markers in LD. Default: 0.05. Min value: value used when generating the tld file.", false, true),
    maxP("p", "maxP", "The maximal p-value to consider. Default: 1e-6.", false, true),
    idColName("id", "idColName", "The name of the variant id column. Default: 'variantId'.", false, true),
    pColName("pn", "pColName", "The name of the p-value column. Default: 'h.intercept.p'.", false, true),
    contigColName("cn", "contigColName", "The name of the contig column. Default: 'contig'.", false, true),
    phenoColName("phn", "phenoColName", "The name of the phenotype column. Ignored if not provided.", false, true),
    separator("s", "separator", "Separator for the columns. Default: '\t'.", false, true),
    out("o", "out", "The file where to write the results.", true, true);

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
    private LdPruningOptions(
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

        for (LdPruningOptions option : values()) {

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
