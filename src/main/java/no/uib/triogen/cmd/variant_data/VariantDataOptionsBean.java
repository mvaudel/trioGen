package no.uib.triogen.cmd.variant_data;

import java.io.File;
import java.util.HashMap;
import java.util.TreeSet;
import no.uib.triogen.io.covariates.SpecificCovariatesFile;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
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

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(VariantDataOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = aLine.getOptionValue(VariantDataOptions.chromosome.opt);

        // The variant ids
        filePath = aLine.getOptionValue(VariantDataOptions.variantId.opt);

        variantFile = new File(filePath);

        if (!variantFile.exists()) {

            throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

        }

        // the trio file
        filePath = aLine.getOptionValue(VariantDataOptions.trio.opt);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = aLine.getOptionValue(VariantDataOptions.phenoFile.opt);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = aLine.getOptionValue(VariantDataOptions.phenoName.opt);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The general covariates
        if (aLine.hasOption(VariantDataOptions.covariate_general.opt)) {

            option = aLine.getOptionValue(VariantDataOptions.covariate_general.opt);

            covariatesGeneral = option.split(",");

        } else {

            covariatesGeneral = new String[0];

        }

        // The specific covariates
        if (aLine.hasOption(VariantDataOptions.covariate_specific.opt)) {

            option = aLine.getOptionValue(VariantDataOptions.covariate_specific.opt);

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
        if (aLine.hasOption(VariantDataOptions.childId.opt)) {

            option = aLine.getOptionValue(VariantDataOptions.childId.opt);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The output file
        filePath = aLine.getOptionValue(VariantDataOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Variant log
        variantLog = aLine.hasOption(VariantDataOptions.variantLog.opt);

    }
}
