package no.uib.triogen.cmd.mendelian_check;

import java.io.File;
import no.uib.triogen.cmd.ld_matrix.LdMatrixOptions;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class MendelianCheckOptionsBean {

    /**
     * The genotypes file.
     */
    public final File genotypesFile;
    /**
     * The chromosome name.
     */
    public final String chromosome;
    /**
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * The trio file.
     */
    public final File trioFile;
    /**
     * The allele frequency threshold.
     */
    public double alleleFrequencyThreshold = 0.005;
    /**
     * File where to write the output.
     */
    public final File destinationFile;
    /**
     * The number of variants to process simultaneously.
     */
    public int nVariants = Runtime.getRuntime().availableProcessors();
    /**
     * The number of days before timeout.
     */
    public int timeOut = 365;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public MendelianCheckOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (MendelianCheckOptions option : MendelianCheckOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = CliUtils.getOptionValue(aLine, MendelianCheckOptions.geno);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = CliUtils.getOptionValue(aLine, LdMatrixOptions.chromosome);

        // The variant ids
        if (CliUtils.hasOption(aLine, MendelianCheckOptions.variantId)) {

            filePath = CliUtils.getOptionValue(aLine, MendelianCheckOptions.variantId);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The trio file
        filePath = CliUtils.getOptionValue(aLine, MendelianCheckOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The allele frequency threshold
        if (CliUtils.hasOption(aLine, MendelianCheckOptions.af)) {
            
            String option = CliUtils.getOptionValue(aLine, MendelianCheckOptions.af);

            try {

                alleleFrequencyThreshold = Double.parseDouble(option);

                if (alleleFrequencyThreshold < 0.0 || alleleFrequencyThreshold > 1.0) {

                    throw new IllegalArgumentException(
                            "Input for allele frequency threshold (" + option + ") must be a number between 0 and 1."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for allele frequency threshold could not be parsed as a number: " + option + "."
                );

            }
        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, MendelianCheckOptions.out);

        destinationFile = new File(filePath);

        File destinationFolder = destinationFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Number of variants to chew in parallel
        if (CliUtils.hasOption(aLine, MendelianCheckOptions.nVariants)) {

            String argString = CliUtils.getOptionValue(aLine, MendelianCheckOptions.nVariants);

            try {

                nVariants = Integer.parseInt(argString);

                if (nVariants <= 0) {

                    throw new IllegalArgumentException(
                            "Input for number of variants must be a strictly positive number."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for number of variants could not be parsed as a number: " + argString + "."
                );

            }
        }

        // Timeout
        if (CliUtils.hasOption(aLine, MendelianCheckOptions.timeOut)) {

            String argString = CliUtils.getOptionValue(aLine, MendelianCheckOptions.timeOut);

            try {

                timeOut = Integer.parseInt(argString);

                if (timeOut <= 0) {

                    throw new IllegalArgumentException(
                            "Input for timeout must be a positive number."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for timeout could not be parsed as a number: " + argString + "."
                );

            }
        }

    }
}
