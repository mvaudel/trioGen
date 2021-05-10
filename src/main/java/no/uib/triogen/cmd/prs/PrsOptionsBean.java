package no.uib.triogen.cmd.prs;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsComputer;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class PrsOptionsBean {

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
     * Name of the effect allele column in the training file.
     */
    public String eaColumn = "tested_allele";
    /**
     * The training file.
     */
    public final File trainingFile;
    /**
     * Name of the variant identifier column in the training file.
     */
    public String snpIdColumn = "variantId";
    /**
     * Name of the phenotype in the pheno file.
     */
    public final String phenoName;
    /**
     * Name of the effect allele column in the training file.
     */
    public Model model = Model.cmf;
    /**
     * Names of the variables.
     */
    public String[] variables = new String[]{"c", "m", "f"};
    /**
     * Pattern for the effect size column.
     */
    public String betaPattern = PrsComputer.MODEL_WILDCARD + ".B" + PrsComputer.VARIABLE_WILDCARD;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * LD R2 threshold.
     */
    public double ldThreshold = 0.05;
    /**
     * Allele frequency threshold.
     */
    public double afThreshold = 0.001;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public PrsOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsOptions option : PrsOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(PrsOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = aLine.getOptionValue(PrsOptions.chromosome.opt);

        // The variant ids
        filePath = aLine.getOptionValue(PrsOptions.trainingFile.opt);

        trainingFile = new File(filePath);

        if (!trainingFile.exists()) {

            throw new IllegalArgumentException("Training file (" + trainingFile + ") not found.");

        }

        // the trio file
        filePath = aLine.getOptionValue(PrsOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = aLine.getOptionValue(PrsOptions.phenoFile.opt);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno name
        phenoName = aLine.getOptionValue(PrsOptions.phenoName.opt);

        // The child id column in the pheno file
        if (aLine.hasOption(PrsOptions.childId.opt)) {

            String option = aLine.getOptionValue(PrsOptions.childId.opt);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // snp id column
        if (aLine.hasOption(PrsOptions.snpId.opt)) {

            snpIdColumn = aLine.getOptionValue(PrsOptions.snpId.opt);

        }

        // ea column
        if (aLine.hasOption(PrsOptions.ea.opt)) {

            eaColumn = aLine.getOptionValue(PrsOptions.ea.opt);

        }

        // The model
        if (aLine.hasOption(PrsOptions.model.opt)) {

            String option = aLine.getOptionValue(PrsOptions.model.opt);

            model = Model.valueOf(option);

        }

        // The variables
        if (aLine.hasOption(PrsOptions.variables.opt)) {

            String option = aLine.getOptionValue(PrsOptions.variables.opt);

            variables = option.split(",");

            if (variables.length != model.betaNames.length) {

                String modelVariables = Arrays.stream(variables).collect(Collectors.joining(","));

                throw new IllegalArgumentException("Found " + variables.length + " variables (" + option + ") where " + model.betaNames.length + " (" + modelVariables + ") expected.");

            }
        }

        // The beta pattern
        if (aLine.hasOption(PrsOptions.betaPattern.opt)) {

            betaPattern = aLine.getOptionValue(PrsOptions.betaPattern.opt);

        }

        // The ld threshold
        if (aLine.hasOption(PrsOptions.ldThreshold.opt)) {

            String option = aLine.getOptionValue(PrsOptions.ldThreshold.opt);

            try {

                ldThreshold = Double.parseDouble(option);

                if (Double.isNaN(ldThreshold) || Double.isFinite(ldThreshold)) {

                    throw new IllegalArgumentException("The value for LD threshold (" + option + ") could not be parsed as a number.");

                }
                if (ldThreshold < 0 || ldThreshold > 1) {

                    throw new IllegalArgumentException("The LD threshold (" + option + ") should be higher than 0 or lower than 1.");

                }
            } catch (Exception e) {

                throw new IllegalArgumentException("The value for LD threshold (" + option + ") could not be parsed as a number.");

            }
        }

        // The af threshold
        if (aLine.hasOption(PrsOptions.afThreshold.opt)) {

            String option = aLine.getOptionValue(PrsOptions.afThreshold.opt);

            try {

                afThreshold = Double.parseDouble(option);

                if (Double.isNaN(afThreshold) || Double.isFinite(afThreshold)) {

                    throw new IllegalArgumentException("The value for allele frequency threshold (" + option + ") could not be parsed as a number.");

                }
                if (afThreshold < 0 || afThreshold > 1) {

                    throw new IllegalArgumentException("The allele frequency threshold (" + option + ") should be higher than 0 or lower than 1.");

                }
            } catch (Exception e) {

                throw new IllegalArgumentException("The value for allele frequency threshold (" + option + ") could not be parsed as a number.");

            }
        }

        // The output file
        filePath = aLine.getOptionValue(PrsOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
