package no.uib.triogen.cmd.ld;

import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesFileType;
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
     * The genotypes file type.
     */
    public GenotypesFileType genotypesFileType = GenotypesFileType.vcf;
    /**
     * the trio file.
     */
    public final File trioFile;
    /**
     * the max distance.
     */
    public int maxDistance = 500000;
    /**
     * The maf threshold.
     */
    public double maf = 0.05;
    /**
     * Boolean indicating whether hard calls should be used.
     */
    public boolean hardCalls = false;
    /**
     * File where to write the output.
     */
    public final String destinationFilePath;
    /**
     * The number of variants to process simultaneously.
     */
    public int nVariants = Runtime.getRuntime().availableProcessors();
    ;
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

        // The genotypes file type
        if (aLine.hasOption(LdMatrixOptions.genoFormat.opt)) {

            String option = aLine.getOptionValue(LdMatrixOptions.genoFormat.opt);
            int genoFormat;

            try {

                genoFormat = Integer.valueOf(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("Genotype file formant could not be parsed. Found: " + option + ". Expected input: " + GenotypesFileType.getCommandLineOptions() + ".");

            }

            genotypesFileType = GenotypesFileType.getGenotypesFileType(genoFormat);

        }

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

        // The maf threshold
        if (aLine.hasOption(LdMatrixOptions.maf.opt)) {
            
            String option = aLine.getOptionValue(LdMatrixOptions.maf.opt);

            try {

                maf = Double.parseDouble(option);

                if (maf < 0.0 || maf > 1.0) {

                    throw new IllegalArgumentException(
                            "Input for maf (" + option + ") must be a number between 0 and 1."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for timeout could not be parsed as a number: " + option + "."
                );

            }
        }

        // Hard calls
        if (aLine.hasOption(LdMatrixOptions.hardCalls.opt)) {

            hardCalls = true;
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

        // Test
        test = aLine.hasOption(LdMatrixOptions.test.opt);

    }
}
