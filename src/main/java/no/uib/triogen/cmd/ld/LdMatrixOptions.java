package no.uib.triogen.cmd.ld;

import java.util.Arrays;
import org.apache.commons.cli.Options;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;

/**
 * Enum of the different options
 *
 * @author Marc Vaudel
 */
public enum LdMatrixOptions {

    geno("g", "geno", "The genotypes file.", true, true),
    genoFormat("gf", "genoFormat", "The genotypes file format to use when iterating the entire file. " + GenotypesFileType.getCommandLineOptions() + ". Default: " + GenotypesFileType.vcf.index + ".", false, true),
    variantId("vi", "variantId", "File listing the variants to include in the analysis. Default: process all variants in the genotypes file.", false, true),
    maxDistance("d", "dist", "The maximum distance in bp to consider around a variant. Default: 500000.", false, true),
    trio("f", "fam", "The trio identifiers file. Can be gzipped or not. Consider including only unrelated samples and controling for admixture.", true, true),
    minR2("r", "minR2", "The minimal ld r2 to report (inclusive). Default: 1e-6.", false, true),
    maf("maf", "mafThreshold", "Minor allele frequency threshold. 0.05 excludes all variants with a maf < 5% among all parents belonging to trios where a phenotype is available for the regression. Default: 0.05.", false, true),
    hardCalls("hc", "hard_calls", "If present, use hard calls instead of dosages.", false, false),
    out("o", "out", "The file where to write the matrix. The extension '.tld' will be added if not present.", true, true),
    nVariants("nv", "nVariants", "The number of variants to process in parallel. Default is the number of cores on the machine.", false, true),
    downstreamLoadingFactor("dlf", "downstreamLoadingFactor", "The downstream loading factor, see documentation for details. Default: 1.1.", false, true),
    upstreamLoadingFactor("ulf", "pustreamLoadingFactor", "The upstream loading factor, see documentation for details. Default: 1.2.", false, true),
    timeOut("z", "timeOut", "The number of days before timeout, default is 365.", false, true),
    test("t", "test", "If present, runs only othe first 1000 variants.", false, false),
    testIteration("ti", "testIteration", "If present, iterate through the variants but do not compute LD calculations.", false, false);

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
    private LdMatrixOptions(
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

        for (LdMatrixOptions option : values()) {

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
