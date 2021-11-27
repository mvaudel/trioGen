package no.uib.triogen.cmd.simple_score;

import java.io.File;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class SimpleScoreOptionsBean {

    /**
     * The genotypes file.
     */
    public final File genotypesFile;
    /**
     * The chromosome name.
     */
    public final String chromosome;
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
    public final File variantFile;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * The number of days before timeout.
     */
    public int timeOut = 365;
    /**
     * Variant log.
     */
    public final boolean variantLog;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public SimpleScoreOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (SimpleScoreOptions option : SimpleScoreOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = CliUtils.getOptionValue(aLine, SimpleScoreOptions.geno);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = CliUtils.getOptionValue(aLine, SimpleScoreOptions.chromosome);

        // The variant ids
        filePath = CliUtils.getOptionValue(aLine, SimpleScoreOptions.variantId);

        variantFile = new File(filePath);

        if (!variantFile.exists()) {

            throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

        }

        // the trio file
        filePath = CliUtils.getOptionValue(aLine, SimpleScoreOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = CliUtils.getOptionValue(aLine, SimpleScoreOptions.phenoFile);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = CliUtils.getOptionValue(aLine, SimpleScoreOptions.phenoName);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The child id column in the pheno file
        if (CliUtils.hasOption(aLine, SimpleScoreOptions.childId)) {

            option = CliUtils.getOptionValue(aLine, SimpleScoreOptions.childId);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, SimpleScoreOptions.out);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Timeout
        if (CliUtils.hasOption(aLine, SimpleScoreOptions.timeOut)) {

            option = CliUtils.getOptionValue(aLine, SimpleScoreOptions.timeOut);

            try {

                timeOut = Integer.parseInt(option);

                if (timeOut <= 0) {

                    throw new IllegalArgumentException(
                            "Input for timeout (" + option + ") must be a positive number."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for timeout could not be parsed as a number: " + option + "."
                );

            }
        }

        // Variant log
        variantLog = CliUtils.hasOption(aLine, SimpleScoreOptions.variantLog);

    }
}
