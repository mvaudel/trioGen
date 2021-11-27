package no.uib.triogen.cmd.association;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.covariates.SpecificCovariatesFile;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.processing.linear_model.LinearModelRunnable;
import no.uib.triogen.utils.cli.CliUtils;
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
    public File variantFile = null;
    /**
     * The max distance to allow between snps.
     */
    public int maxDistance = 500000;
    /**
     * The allele frequency threshold.
     */
    public double alleleFrequencyThreshold = 0.005;
    /**
     * List of the names of the models to use.
     */
    public String[] modelNames = Model.getDefaultOption();
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
     * Variant log.
     */
    public final boolean variantLog;

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

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = CliUtils.getOptionValue(aLine, LinearModelOptions.geno);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Genotypes file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = CliUtils.getOptionValue(aLine, LinearModelOptions.chromosome);

        // The variant ids
        if (CliUtils.hasOption(aLine, LinearModelOptions.variantId)) {

            filePath = CliUtils.getOptionValue(aLine, LinearModelOptions.variantId);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The max distance
        if (CliUtils.hasOption(aLine, LinearModelOptions.maxDistance)) {

            String stringValue = CliUtils.getOptionValue(aLine, LinearModelOptions.maxDistance);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance < 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a positive integer.");

            }
        }

        // The allele frequency threshold
        if (CliUtils.hasOption(aLine, LinearModelOptions.afThreshold)) {

            String option = CliUtils.getOptionValue(aLine, LinearModelOptions.afThreshold);

            try {

                alleleFrequencyThreshold = Double.parseDouble(option);

                if (alleleFrequencyThreshold < 0.0 || alleleFrequencyThreshold > 1.0) {

                    throw new IllegalArgumentException(
                            "Input for allele frequency (" + option + ") must be a number between 0 and 1."
                    );
                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for timeout could not be parsed as a number: " + option + "."
                );

            }
        }

        // the trio file
        filePath = CliUtils.getOptionValue(aLine, LinearModelOptions.trio);

        trioFile = new File(filePath);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file (" + trioFile + ") not found.");

        }

        // The pheno file
        filePath = CliUtils.getOptionValue(aLine, LinearModelOptions.phenoFile);

        phenotypesFile = new File(filePath);

        if (!phenotypesFile.exists()) {

            throw new IllegalArgumentException("Phenotype file (" + phenotypesFile + ") not found.");

        }

        // The pheno columns
        String option = CliUtils.getOptionValue(aLine, LinearModelOptions.phenoName);

        phenoNames = option.split(",");

        if (option.length() == 0 || phenoNames.length == 0) {

            throw new IllegalArgumentException("No phenotype name found.");

        }

        // The general covariates
        if (CliUtils.hasOption(aLine, LinearModelOptions.covariate_general)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.covariate_general);

            covariatesGeneral = option.split(",");

        } else {

            covariatesGeneral = new String[0];

        }

        // The specific covariates
        if (CliUtils.hasOption(aLine, LinearModelOptions.covariate_specific)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.covariate_specific);

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
        if (CliUtils.hasOption(aLine, LinearModelOptions.childId)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.childId);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The models
        if (CliUtils.hasOption(aLine, LinearModelOptions.model)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.model);

            modelNames = option.split(",");

            if (option.length() == 0 || modelNames.length == 0) {

                throw new IllegalArgumentException("No model found.");

            }

            HashSet<String> implementedModels = Arrays.stream(Model.values())
                    .map(
                            model -> model.name()
                    )
                    .collect(
                            Collectors.toCollection(HashSet::new)
                    );
            String missingModels = Arrays.stream(modelNames)
                    .filter(
                            model -> !implementedModels.contains(model)
                    )
                    .collect(
                            Collectors.joining(",")
                    );

            if (missingModels.length() > 0) {

                throw new IllegalArgumentException("Models not implemented: " + missingModels + ". Please use one of the following models: " + Model.getCommandLineOptions() + ".");

            }
        }

        HashSet<String> includedModels = Arrays.stream(modelNames)
                .collect(
                        Collectors.toCollection(HashSet::new)
                );
        Arrays.stream(Model.values())
                .forEach(
                        model -> model.setParentModels(includedModels)
                );

        // Inclusion of cases where no regression can be done
        if (CliUtils.hasOption(aLine, LinearModelOptions.x0)) {

            LinearModelRunnable.x0 = true;

        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, LinearModelOptions.out);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

        // Number of variants to chew in parallel
        if (CliUtils.hasOption(aLine, LinearModelOptions.nVariants)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.nVariants);

            try {

                nVariants = Integer.parseInt(option);

                if (timeOut <= 0) {

                    throw new IllegalArgumentException(
                            "Input for number of variants (" + option + ") must be a strictly positive number."
                    );

                }

            } catch (Exception e) {

                e.printStackTrace();

                throw new IllegalArgumentException(
                        "Input for number of variants could not be parsed as a number: " + option + "."
                );

            }
        }

        // Timeout
        if (CliUtils.hasOption(aLine, LinearModelOptions.timeOut)) {

            option = CliUtils.getOptionValue(aLine, LinearModelOptions.timeOut);

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
        variantLog = CliUtils.hasOption(aLine, LinearModelOptions.variantLog);

    }
}
