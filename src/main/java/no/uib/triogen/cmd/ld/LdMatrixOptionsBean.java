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
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * the trio file.
     */
    public final File trioFile;
    /**
     * the max distance.
     */
    public int maxDistance = 1000000;
    /**
     * The minimal ld r2 to report (inclusive).
     */
    public double minR2 = 1e-6;
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
    /**
     * The downstream loading factor.
     */
    public double downstreamLoadingFactor = 1.05;
    /**
     * The upstream loading factor.
     */
    public double upstreamLoadingFactor = 1.5;
    /**
     * The number of days before timeout.
     */
    public int timeOut = 365;
    /**
     * Test mode.
     */
    public final boolean test;
    /**
     * Iteration test mode.
     */
    public final boolean testIteration;

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

        // The variant ids
        if (aLine.hasOption(LdMatrixOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(LdMatrixOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
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
                        "Input for maf could not be parsed as a number: " + option + "."
                );

            }
        }

        // The minimal R2
        if (aLine.hasOption(LdMatrixOptions.minR2.opt)) {
            
            String option = aLine.getOptionValue(LdMatrixOptions.minR2.opt);

            try {

                minR2 = Double.parseDouble(option);

                if (minR2 <= 0.0 || minR2 >= 1.0) {

                    throw new IllegalArgumentException(
                            "Input for maf (" + option + ") must be a number between 0 and 1 (both excluded)."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for minR2 could not be parsed as a number: " + option + "."
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

        // The downstream loading factor
        if (aLine.hasOption(LdMatrixOptions.downstreamLoadingFactor.opt)) {

            String argString = aLine.getOptionValue(LdMatrixOptions.downstreamLoadingFactor.opt);

            try {

                downstreamLoadingFactor = Double.parseDouble(argString);

                if (downstreamLoadingFactor < 1) {

                    throw new IllegalArgumentException(
                            "Input for downstream loading factor must be higher or equal to one."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for downstream loading factor could not be parsed as a number: " + argString + "."
                );
            }
        }

        // The upstream loading factor
        if (aLine.hasOption(LdMatrixOptions.upstreamLoadingFactor.opt)) {

            String argString = aLine.getOptionValue(LdMatrixOptions.upstreamLoadingFactor.opt);

            try {

                upstreamLoadingFactor = Double.parseDouble(argString);

                if (upstreamLoadingFactor < 1) {

                    throw new IllegalArgumentException(
                            "Input for upstream loading factor must be higher or equal to one."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for upstream loading factor could not be parsed as a number: " + argString + "."
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

        // Test iteration
        testIteration = aLine.hasOption(LdMatrixOptions.testIteration.opt);

    }
}
