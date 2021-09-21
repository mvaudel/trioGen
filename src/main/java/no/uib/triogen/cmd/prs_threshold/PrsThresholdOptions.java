package no.uib.triogen.cmd.prs_threshold;

import java.util.Arrays;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import org.apache.commons.cli.Options;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsPruner;
import static no.uib.triogen.utils.Utils.CHROMOSOME_WILDCARD;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum PrsThresholdOptions {

    geno("g", "geno", "The genotypes file. Wildcard: " + CHROMOSOME_WILDCARD + " for variable name.", true, true),
    phenoFile("p", "phenoFile", "The phenotypes file.", true, true),
    phenoName("pn", "phenoName", "The name of the phenotype.", true, true),
    variantId("vi", "variantId", "File listing the variants to include in the analysis. Default: process all variants in the genotypes file.", false, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not.", true, true),
    scoreFile("s", "scoreFile", "Score details as obtained from the PrsTrain command.", true, true),
    model("m", "model", "Names of the model to use. Default: " + Model.cmf.name() + ". Available: " + Model.getCommandLineOptions() + ".", false, true),
    variables("v", "variables", "Names of the variables to use, need to be in the same order as specified in the model. Default: c,m,f.", false, true),
    nBins("nb", "nBins", "The number of bins to use when comparing scores and phenotyes. Ignored if zero. Default: 10.", false, true),
    betaPattern("bp", "betaPattern", "Pattern for the effect size column. Wildcard: " + PrsPruner.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD + "'.", false, true),
    sePattern("sp", "sePattern", "Pattern for the standard error column. Wildcard: " + PrsPruner.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD + ".se'.", false, true),
    afThreshold("af", "afThreshold", "Allele frequency threshold. 0.005 excludes all alleles of variants with frequency < 0.5% or > 99.5%. Default: no threshold.", false, true),
    afColumn("ac", "afColumn", "Allele frequency clumn. Required if an allele frequency threshold is set.", false, true),
    out("o", "out", "The file where to write the results.", true, true),
    bgenIndexFolder("bif", "bgenIndexFolder", "The folder for the bgen index files. Default: the folder where the bgen files are located.", false, true);

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
    private PrsThresholdOptions(
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

        for (PrsThresholdOptions option : values()) {

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
