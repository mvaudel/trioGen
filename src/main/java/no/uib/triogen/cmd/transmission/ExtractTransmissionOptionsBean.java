package no.uib.triogen.cmd.transmission;

import java.io.File;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class ExtractTransmissionOptionsBean {

    /**
     * the vcf file.
     */
    public final File vcfFile;
    /**
     * the trio file.
     */
    public final File trioFile;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * The number of threads to use.
     */
    public int nThreads = 1;
    /**
     * The number of days before timeout.
     */
    public int timeOut = 365;
    /**
     * Test mode.
     */
    public final boolean test;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public ExtractTransmissionOptionsBean(CommandLine aLine) {

        // Check that mandatory options are provided
        for (ExtractTransmissionOptions option : ExtractTransmissionOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // Vcf
        String filePath = aLine.getOptionValue(ExtractTransmissionOptions.vcf.opt);

        vcfFile = new File(filePath);

        if (!vcfFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Vcf file (" + vcfFile.getParent() + ") not found.");

        }

        // trio
        filePath = aLine.getOptionValue(ExtractTransmissionOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile.getParent() + ") not found.");

        }

        // Output
        filePath = aLine.getOptionValue(ExtractTransmissionOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output file (" + destinationFile.getParent() + ") not found.");

        }

        // Number of threads
        if (aLine.hasOption(ExtractTransmissionOptions.threads.opt)) {

            String argString = aLine.getOptionValue(ExtractTransmissionOptions.threads.opt);

            try {

                nThreads = Integer.parseInt(argString);

                if (nThreads <= 0) {

                    throw new IllegalArgumentException(
                            "Input for number of threads must be a positive number."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for number of threads could not be parsed as a number: " + argString + "."
                );

            }
        }

        // Timeout
        if (aLine.hasOption(ExtractTransmissionOptions.timeOut.opt)) {

            String argString = aLine.getOptionValue(ExtractTransmissionOptions.timeOut.opt);

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

        // Test
        test = aLine.hasOption(ExtractTransmissionOptions.test.opt);

    }
}
