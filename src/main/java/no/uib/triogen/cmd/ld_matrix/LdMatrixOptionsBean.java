package no.uib.triogen.cmd.ld_matrix;

import java.io.File;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LdMatrixOptionsBean {

    /**
     * The genotypes file.
     */
    public final File genotypesFile;
    /**
     * The chromosome name.
     */
    public final String chromosome;
    /**
     * The max distance between the snp and the target snp.
     */
    public int maxDistance = 500000;
    /**
     * The trio file.
     */
    public final File trioFile;
    /**
     * The minimal ld r2 to report (inclusive).
     */
    public double minR2 = 1e-6;
    /**
     * The allele frequency threshold.
     */
    public double alleleFrequencyThreshold = 0.001;
    /**
     * File where to write the output.
     */
    public final String destinationFilePath;
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
    public LdMatrixOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (LdMatrixOptions option : LdMatrixOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = CliUtils.getOptionValue(aLine, LdMatrixOptions.geno);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = CliUtils.getOptionValue(aLine, LdMatrixOptions.chromosome);

        // The trio file
        filePath = CliUtils.getOptionValue(aLine, LdMatrixOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The max distance
        if (CliUtils.hasOption(aLine, LdMatrixOptions.maxDistance)) {

            String stringValue = CliUtils.getOptionValue(aLine, LdMatrixOptions.maxDistance);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance <= 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a stricly positive integer.");

            }
        }

        // The minimal R2
        if (CliUtils.hasOption(aLine, LdMatrixOptions.minR2)) {
            
            String option = CliUtils.getOptionValue(aLine, LdMatrixOptions.minR2);

            try {

                minR2 = Double.parseDouble(option);

                if (minR2 <= 0.0 || minR2 >= 1.0) {

                    throw new IllegalArgumentException(
                            "Input for minimal R2 (" + option + ") must be a number between 0 and 1 (both excluded)."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for minR2 could not be parsed as a number: " + option + "."
                );

            }
        }

        // The allele frequency threshold
        if (CliUtils.hasOption(aLine, LdMatrixOptions.af)) {

            String option = CliUtils.getOptionValue(aLine, LdMatrixOptions.af);

            try {

                alleleFrequencyThreshold = Double.parseDouble(option);

                if (alleleFrequencyThreshold < 0.0 || alleleFrequencyThreshold > 1.0) {

                    throw new IllegalArgumentException(
                            "Input for allele frequency (" + option + ") must be a number between 0 and 1."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for timeout could not be parsed as a number: " + option + "."
                );

            }
        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, LdMatrixOptions.out);
        
        if (filePath.endsWith(".tld")) {
            
            filePath = filePath.substring(0, filePath.length() - 4);
            
        }

        destinationFilePath = filePath;

        File destinationFolder = (new File(destinationFilePath)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Number of variants to chew in parallel
        if (CliUtils.hasOption(aLine, LdMatrixOptions.nVariants)) {

            String argString = CliUtils.getOptionValue(aLine, LdMatrixOptions.nVariants);

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
        if (CliUtils.hasOption(aLine, LdMatrixOptions.timeOut)) {

            String argString = CliUtils.getOptionValue(aLine, LdMatrixOptions.timeOut);

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
