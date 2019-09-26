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
     * Stem of the file where to write the output.
     */
    public final String destinationStem;
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
    public ExtractTransmissionOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (ExtractTransmissionOptions option : ExtractTransmissionOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // Vcf
        String filePath = aLine.getOptionValue(ExtractTransmissionOptions.vcf.opt);

        vcfFile = new File(filePath);

        if (!vcfFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + vcfFile + ") not found.");

        }

        // trio
        filePath = aLine.getOptionValue(ExtractTransmissionOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // Output
        filePath = aLine.getOptionValue(ExtractTransmissionOptions.out.opt);

        destinationStem = filePath;
        
        File destinationFolder = (new File(destinationStem)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

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
