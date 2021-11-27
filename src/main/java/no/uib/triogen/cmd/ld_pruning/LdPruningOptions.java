package no.uib.triogen.cmd.ld_pruning;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.utils.Utils;
import no.uib.triogen.utils.cli.CliOption;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum LdPruningOptions implements CliOption {

    results("res", "results", "The results file to prune.", true, true),
    ldMatrix("l", "ld_matrix", "The ld matrix file as generated using the LdMatrix command. If LD matrix files are computed per contig, replace the chromosome name with '" + Utils.CHROMOSOME_WILDCARD + "'. Ignored if not provided.", false, true),
    build("b", "build", "The build to use when querying Ensembl as a number 37: grch37, 38: grch38. Default: 37.", false, true),
    ensemblPopulation("ep", "ensembl_population", "The reference population for Ensembl. See https://rest.ensembl.org/documentation/info/variation_populations for details. Ensembl is not used if not provided.", false, true),
    ldlinkPopulation("lp", "ldlink_population", "The reference population to use for LDlink. LDlink is not used if not provided.", false, true),
    ldlinkToken("lt", "ldlink_token", "The token to use for LDlink. Mandatory when using LDlink.", false, true),
    minR2("r", "min_r2", "The minimal ld r2 to consider two markers in LD. Default: 0.05.", false, true),
    maxP("p", "max_p", "The maximal p-value to consider. Default: 1e-6.", false, true),
    idColName("id", "id_col_name", "The name of the variant id column. Default: 'variantId'.", false, true),
    rsidColName("rs", "rsid_col_name", "The name of the rsid column. Mandatory when using Ensembl or LDlink", false, true),
    pColName("pn", "p_col_name", "The name of the p-value column. Default: 'h.intercept.p'.", false, true),
    contigColName("cn", "contig_col_name", "The name of the contig column. Default: 'contig'.", false, true),
    phenoColName("phn", "pheno_col_name", "The name of the phenotype column. Ignored if not provided.", false, true),
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

    @Override
    public String getOption() {
        
        return opt;
        
    }

    @Override
    public String getLongOption() {
        
        return longOpt;
        
    }
}
