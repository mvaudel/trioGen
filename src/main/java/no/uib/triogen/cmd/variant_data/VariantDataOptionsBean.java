package no.uib.triogen.cmd.variant_data;

import java.io.File;
import java.util.HashMap;
import java.util.TreeSet;
import no.uib.triogen.io.covariates.SpecificCovariatesFile;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class VariantDataOptionsBean {

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
     * Names of the covariates to use for all the phenotypes.
     */
    public final String[] covariatesGeneral;
    /**
     * Map of the covariates to use for specific phenotypes.
     */
    public final HashMap<String, TreeSet<String>> covariatesSpecific;
    /**
     * The file listing the variants to process.
     */
    public final File variantFile;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
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
    public VariantDataOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (VariantDataOptions option : VariantDataOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = CliUtils.getOptionValue(aLine, VariantDataOptions.geno);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = CliUtils.getOptionValue(aLine, VariantDataOptions.chromosome);

        // The variant ids
        filePath = CliUtils.getOptionValue(aLine, VariantDataOptions.variantId);

        variantFile = new File(filePath);

        if (!variantFile.exists()) {

            throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

        }

        // the trio file
        filePath = CliUtils.getOptionValue(aLine, VariantDataOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = CliUtils.getOptionValue(aLine, VariantDataOptions.phenoFile);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = CliUtils.getOptionValue(aLine, VariantDataOptions.phenoName);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The general covariates
        if (CliUtils.hasOption(aLine, VariantDataOptions.covariate_general)) {

            option = CliUtils.getOptionValue(aLine, VariantDataOptions.covariate_general);

            covariatesGeneral = option.split(",");

        } else {

            covariatesGeneral = new String[0];

        }

        // The specific covariates
        if (CliUtils.hasOption(aLine, VariantDataOptions.covariate_specific)) {

            option = CliUtils.getOptionValue(aLine, VariantDataOptions.covariate_specific);

            File covariatesFile = new File(option);

            if (!covariatesFile.exists()) {

                throw new IllegalArgumentException("Covariates file " + option + " not found.");

            }

            covariatesSpecific = SpecificCovariatesFile.praseCovariates(covariatesFile);

        } else {

            covariatesSpecific = new HashMap<>(0);

        }

        for (String phenoName : phenoNames) {

            if (!covariatesSpecific.containsKey(phenoName)) {

                covariatesSpecific.put(phenoName, new TreeSet<>());

            }
        }

        // The child id column in the pheno file
        if (CliUtils.hasOption(aLine, VariantDataOptions.childId)) {

            option = CliUtils.getOptionValue(aLine, VariantDataOptions.childId);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, VariantDataOptions.out);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Variant log
        variantLog = CliUtils.hasOption(aLine, VariantDataOptions.variantLog);

    }
}
