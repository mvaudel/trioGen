package no.uib.triogen.cmd.summary;

import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class SummaryOptionsBean {

    /**
     * The genotypes file.
     */
    public final File genotypesFile;
    /**
     * The genotypes file type.
     */
    public GenotypesFileType genotypesFileType = GenotypesFileType.sangerVCF;
    /**
     * the trio file.
     */
    public final File trioFile;
    /**
     * Stem of the file where to write the output.
     */
    public final String destinationStem;
    /**
     * The number of variants to process simultaneously.
     */
    public int nVariants = 8;
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
    public SummaryOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (SummaryOptions option : SummaryOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(SummaryOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The genotypes file type
        if (aLine.hasOption(SummaryOptions.genoFormat.opt)) {

            String option = aLine.getOptionValue(SummaryOptions.genoFormat.opt);
            int genoFormat;

            try {

                genoFormat = Integer.valueOf(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("Genotype file formant could not be parsed. Found: " + option + ". Expected input: " + GenotypesFileType.getCommandLineOptions() + ".");

            }

            genotypesFileType = GenotypesFileType.getGenotypesFileType(genoFormat);

        }

        // the trio file
        filePath = aLine.getOptionValue(SummaryOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The output file
        filePath = aLine.getOptionValue(SummaryOptions.out.opt);

        destinationStem = filePath;

        File destinationFolder = (new File(destinationStem)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Number of variants to chew in parallel
        if (aLine.hasOption(SummaryOptions.nVariants.opt)) {

            String argString = aLine.getOptionValue(SummaryOptions.nVariants.opt);

            try {

                nVariants = Integer.parseInt(argString);

                if (timeOut <= 0) {

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
        if (aLine.hasOption(SummaryOptions.timeOut.opt)) {

            String argString = aLine.getOptionValue(SummaryOptions.timeOut.opt);

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
        test = aLine.hasOption(SummaryOptions.test.opt);

    }
}