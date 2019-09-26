package no.uib.triogen.cmd.transmission;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.Utils.lineSeparator;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum ExtractTransmissionOptions {

    vcf("g", "geno", "The vcf file. Can be gzipped or not.", true, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    out("o", "out", "stem of the files where to write the scores.", true, true),
    timeOut("z", "timeOut", "Number of days before timeout, default is 365.", false, true),
    test("t", "test", "If present, runs only one variant in one trio.", false, false);

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
    private ExtractTransmissionOptions(String opt, String longOpt, String description, boolean mandatory, boolean hasArg) {
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
    public static void createOptionsCLI(Options options) {

        for (ExtractTransmissionOptions option : values()) {

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
