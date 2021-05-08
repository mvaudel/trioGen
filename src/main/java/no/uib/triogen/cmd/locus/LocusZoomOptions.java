package no.uib.triogen.cmd.locus;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum LocusZoomOptions {

    ldMatrix("ld", "ldMatrix", "The file containing the ldMatrix. 'tld' file generated using the 'ldMatrix' command.", true, true),
    results("res", "results", "The results file of the LinearModel command,", true, true),
    variantId("vi", "variantId", "File listing the variants to query.", true, true),
    phenoName("pn", "phenoName", "The name of the phenotype to export the locus zoom data on. Example: pheno1. Default: all phenotypes found for the given variants.", false, true),
    maxDistance("d", "dist", "The maximum distance in bp to consider around a variant. Should be below the distance used to create the ld matrix and, if gene mapping is conducted, below 2.5e6. Default: 1000000.", false, true),
    buildNumber("b", "buildNumber", "The number of the build to use if gene mapping is conducted. e.g. 38 for GRCh38. Default: 38.", false, true),
    out("o", "out", "The file where to write the data needed to build a locus zoom plot.", true, true),
    geneCoordinates("g", "geneCoordinates", "The file where to write the gene coordinates. If none provided gene mapping will be skipped.", false, true),
    log("log", "log", "The file where to write the log. Default: next to the output.", false, true);

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
    private LocusZoomOptions(
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

        for (LocusZoomOptions option : values()) {

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
