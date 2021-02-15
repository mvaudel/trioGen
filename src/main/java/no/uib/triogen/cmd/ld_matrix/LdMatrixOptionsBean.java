package no.uib.triogen.cmd.ld_matrix;

import java.io.File;
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

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(LdMatrixOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = aLine.getOptionValue(LdMatrixOptions.chromosome.opt);

        // The trio file
        filePath = aLine.getOptionValue(LdMatrixOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The max distance
        if (aLine.hasOption(LdMatrixOptions.maxDistance.opt)) {

            String stringValue = aLine.getOptionValue(LdMatrixOptions.maxDistance.opt);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance <= 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a stricly positive integer.");

            }
        }

        // The minimal R2
        if (aLine.hasOption(LdMatrixOptions.minR2.opt)) {
            
            String option = aLine.getOptionValue(LdMatrixOptions.minR2.opt);

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

        // The output file
        filePath = aLine.getOptionValue(LdMatrixOptions.out.opt);
        
        if (filePath.endsWith(".tld")) {
            
            filePath = filePath.substring(0, filePath.length() - 4);
            
        }

        destinationFilePath = filePath;

        File destinationFolder = (new File(destinationFilePath)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Number of variants to chew in parallel
        if (aLine.hasOption(LdMatrixOptions.nVariants.opt)) {

            String argString = aLine.getOptionValue(LdMatrixOptions.nVariants.opt);

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
        if (aLine.hasOption(LdMatrixOptions.timeOut.opt)) {

            String argString = aLine.getOptionValue(LdMatrixOptions.timeOut.opt);

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
