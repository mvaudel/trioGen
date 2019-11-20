package no.uib.triogen.cmd.association;

import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import no.uib.triogen.processing.association.linear_model.LinearModelRunnable;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LinearModelOptionsBean {

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
     * The phenotypes file.
     */
    public final File phenotypesFile;
    /**
     * Names of the phenotypes to use for association.
     */
    public final String[] phenoNames;
    /**
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
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
    public LinearModelOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (LinearModelOptions option : LinearModelOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(LinearModelOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The genotypes file type
        if (aLine.hasOption(LinearModelOptions.genoFormat.opt)) {

            String option = aLine.getOptionValue(LinearModelOptions.genoFormat.opt);
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
        if (aLine.hasOption(LinearModelOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(LinearModelOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // the trio file
        filePath = aLine.getOptionValue(LinearModelOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = aLine.getOptionValue(LinearModelOptions.phenoFile.opt);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = aLine.getOptionValue(LinearModelOptions.phenoName.opt);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The child id column in the pheno file
        if (aLine.hasOption(LinearModelOptions.childId.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.childId.opt);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // Inclusion of cases where no regression can be done
        if (aLine.hasOption(LinearModelOptions.x0.opt)) {

                LinearModelRunnable.x0 = true;

        }

        // The output file
        filePath = aLine.getOptionValue(LinearModelOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Number of variants to chew in parallel
        if (aLine.hasOption(LinearModelOptions.nVariants.opt)) {

            String argString = aLine.getOptionValue(LinearModelOptions.nVariants.opt);

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
        if (aLine.hasOption(LinearModelOptions.timeOut.opt)) {

            String argString = aLine.getOptionValue(LinearModelOptions.timeOut.opt);

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
        test = aLine.hasOption(LinearModelOptions.test.opt);

    }
}
