package no.uib.triogen.cmd.prs_prune;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsPruner;
import no.uib.triogen.utils.Utils;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum PrsPruneOptions {

    trainingFile("t", "trainingFile", "File listing the summary statistics to use for training.", true, true),
    ldMatrix("l", "ldMatrix", "The ld matrix file as generated using the LdMatrix command. If LD matrix files are computed per contig, replace the chromosome name with '" + Utils.CHROMOSOME_WILDCARD + "'.", true, true),
    snpId("sc", "snpId", "Name of the variant identifier column in the training file. Default: 'variantId'.", false, true),
    chrColumn("cc", "chrColumn", "Name of the contig column in the training file. Default: 'contig'.", true, true),
    posColumn("pc", "posColumn", "Name of the position column in the training file. Default: 'position'.", true, true),
    refColumn("rc", "refColumn", "Name of the reference allele column in the training file. Default: 'otherAllele'.", true, true),
    eaColumn("ec", "eaColumn", "Name of the effect allele column in the training file. Default: 'testedAllele'.", true, true),
    variables("v", "variables", "Names of the variables to use, need to be in the same order as specified in the model. Default: c,m,f.", false, true),
    betaPattern("bp", "betaPattern", "Pattern for the effect size column. Wildcard: " + PrsPruner.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD + "'.", false, true),
    sePattern("sp", "sePattern", "Pattern for the standard error column. Wildcard: " + PrsPruner.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD + ".se'.", false, true),
    pPattern("pp", "pPattern", "Pattern for the p-value column. Wildcard: " + PrsPruner.VARIABLE_WILDCARD + " for variable name. Default: '" + Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD + "'.p.", false, true),
    out("o", "out", "The file where to write the results.", true, true),
    ldLocusThreshold("ldl", "ldLocusThreshold", "LD R2 value over which two hits cannot be considered independent. Default: '0.05'.", false, true),
    ldTopHitThreshold("ldh", "ldTopHitThreshold", "LD R2 value over which two hits are considered identical. Default: '0.8'.", false, true),
    pValueThreshold("pv", "pValueThreshold", "The highest p-value to consider. Default: '0.05'.", false, true),
    nSnpPerLocusThreshold("ldn", "nSnpPerLocusThreshold", "The minimal number of variants required for a locus. Default: '5'.", false, true);

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
    private PrsPruneOptions(
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

        for (PrsPruneOptions option : values()) {

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
