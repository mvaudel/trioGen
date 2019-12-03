package no.uib.triogen.cmd.results;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.lineSeparator;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum ExtractOptions {

    input("i", "input", "The results file.", true, true),
    category("cat", "category", "The categories columns as comma separated list, one file will be produced per level. Example: pheno,variantId. Default: no category.", false, true),
    value("val", "value", "The columns to include in the results file as comma separated list. Example: h_B1,h_B1_se,h_B1_p. Default: all columns.", false, true),
    output("o", "output", "Stem where to write the output", true, true);

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
    private ExtractOptions(
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

        for (ExtractOptions option : values()) {

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
        output.append(lineSeparator)
                .append(lineSeparator);
        
        output.append("-").append(String.format(formatter, "h (--help)")).append(" ").append("Shows a brief help message.").append(lineSeparator);
        output.append("-").append(String.format(formatter, "v (--version)")).append(" ").append("Shows the version of the tool.").append(lineSeparator);

        output.append(lineSeparator)
                .append(lineSeparator);
        output.append("Mandatory Options:");
        output.append(lineSeparator)
                .append(lineSeparator);

        Arrays.stream(values())
                .filter(option -> option.mandatory)
                .forEach(option -> output.append("-").append(String.format(formatter, option.opt + " (--" + option.longOpt + ")")).append(" ").append(option.description).append(lineSeparator));

        output.append(lineSeparator)
                .append(lineSeparator);
        output.append("Additional Options:");
        output.append(lineSeparator)
                .append(lineSeparator);

        Arrays.stream(values())
                .filter(option -> !option.mandatory)
                .forEach(option -> output.append("-").append(String.format(formatter, option.opt + " (--" + option.longOpt + ")")).append(" ").append(option.description).append(lineSeparator));

        return output.toString();
    }
}
