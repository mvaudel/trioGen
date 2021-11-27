package no.uib.triogen.cmd.simple_score;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.utils.cli.CliOption;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum SimpleScoreOptions implements CliOption {

    geno("g", "geno", "The genotypes file.", true, true),
    chromosome("c", "chromosome", "The chromosome name.", true, true),
    variantId("vi", "variantId", "File listing the variants to include in the analysis along with the weights.", true, true),
    phenoFile("p", "phenoFile", "The phenotypes file.", true, true),
    childId("id", "childId", "The name of the column containing the child id. Default: child_SentrixID.", false, true),
    phenoName("pn", "phenoName", "List of the names of the phenotypes to include in the export. Example: pheno1,pheno2.", true, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    out("o", "out", "The file where to write the results.", true, true),
    timeOut("z", "timeOut", "The number of days before timeout, default is 365.", false, true),
    variantLog("vl", "variantLog", "If present, writes a log for every variant next to the results file.", false, false);

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
    private SimpleScoreOptions(
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

        for (SimpleScoreOptions option : values()) {

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

    @Override
    public String getOption() {
        
        return opt;
        
    }

    @Override
    public String getLongOption() {
        
        return longOpt;
        
    }
}
