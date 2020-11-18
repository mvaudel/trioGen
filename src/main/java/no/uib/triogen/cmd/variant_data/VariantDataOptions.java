package no.uib.triogen.cmd.variant_data;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.model.trio_genotypes.Model;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum VariantDataOptions {

    geno("g", "geno", "The genotypes file.", true, true),
    genoFormat("gf", "genoFormat", "The genotypes file format to use when iterating the entire file. " + GenotypesFileType.getCommandLineOptions() + ". Default: " + GenotypesFileType.vcf.index + ".", false, true),
    variantId("vi", "variantId", "File listing the variants to include in the analysis.", true, true),
    phenoFile("p", "phenoFile", "The phenotypes file.", true, true),
    childId("id", "childId", "The name of the column containing the child id. Default: child_SentrixID.", false, true),
    phenoName("pn", "phenoName", "List of the names of the phenotypes to include in the analysis. Example: pheno1,pheno2.", true, true),
    covariate_general("cg", "covariate_general", "List of the names of the covariates to use for all phenotypes. Example: pc1,pc2,pc3,pc4,pc5,pc6.", false, true),
    covariate_specific("cs", "covariate_specific", "File containing the names of the covariates to use for a specific phenotype in json format.", false, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    out("o", "out", "The file where to write the results.", true, true),
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
    private VariantDataOptions(
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

        for (VariantDataOptions option : values()) {

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
