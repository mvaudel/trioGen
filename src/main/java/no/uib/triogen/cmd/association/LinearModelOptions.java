package no.uib.triogen.cmd.association;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.Utils.lineSeparator;
import no.uib.triogen.io.genotypes.GenotypesFileType;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum LinearModelOptions {

    geno("g", "geno", "The genotypes file.", true, true),
    genoFormat("gf", "genoFormat", "The genotypes file format. " + GenotypesFileType.getCommandLineOptions(), true, true),
    phenoFile("p", "phenoFile", "The phenotypes file.", true, true),
    childId("id", "childId", "The name of the column containing the child id. Default: child_SentrixID.", false, true),
    phenoName("pn", "phenoName", "List of the names of the phenotypes in the phenotype file. Example: pheno1,pheno2.", true, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    x0("x0", "x0", "If present the association results will only be reported when multiple values of x are available for the regression.", false, false),
    out("o", "out", "The file where to write the results.", true, true),
    nVariants("nv", "nVariants", "The number of variants to process in parallel. Default is 8.", false, true),
    timeOut("z", "timeOut", "The number of days before timeout, default is 365.", false, true),
    test("t", "test", "If present, runs only othe first 1000 variants.", false, false);

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
    private LinearModelOptions(
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

        for (LinearModelOptions option : values()) {

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
