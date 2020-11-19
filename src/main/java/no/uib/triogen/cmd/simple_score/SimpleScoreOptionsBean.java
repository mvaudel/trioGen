package no.uib.triogen.cmd.simple_score;

import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
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
     * The genotypes file type.
     */
    public final GenotypesFileType genotypesFileType = GenotypesFileType.vcf;
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

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(SimpleScoreOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The variant ids

            filePath = aLine.getOptionValue(SimpleScoreOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }

        // the trio file
        filePath = aLine.getOptionValue(SimpleScoreOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = aLine.getOptionValue(SimpleScoreOptions.phenoFile.opt);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = aLine.getOptionValue(SimpleScoreOptions.phenoName.opt);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The child id column in the pheno file
        if (aLine.hasOption(SimpleScoreOptions.childId.opt)) {

            option = aLine.getOptionValue(SimpleScoreOptions.childId.opt);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The output file
        filePath = aLine.getOptionValue(SimpleScoreOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Timeout
        if (aLine.hasOption(SimpleScoreOptions.timeOut.opt)) {

            option = aLine.getOptionValue(SimpleScoreOptions.timeOut.opt);

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
        variantLog = aLine.hasOption(SimpleScoreOptions.variantLog.opt);

    }
}
