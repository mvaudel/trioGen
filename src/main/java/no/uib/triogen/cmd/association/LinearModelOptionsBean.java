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

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotypes file
        String filePath = aLine.getOptionValue(LinearModelOptions.geno.opt);

        genotypesFile = new File(filePath);

        if (!genotypesFile.exists()) {

            throw new IllegalArgumentException("Genotypes file (" + genotypesFile + ") not found.");

        }

        // The chromosome name
        chromosome = aLine.getOptionValue(LinearModelOptions.chromosome.opt);

        // The variant ids
        if (aLine.hasOption(LinearModelOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(LinearModelOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The max distance
        if (aLine.hasOption(LinearModelOptions.maxDistance.opt)) {

            String stringValue = aLine.getOptionValue(LinearModelOptions.maxDistance.opt);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance < 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a positive integer.");

            }
        }

        // The allele frequency threshold
        if (aLine.hasOption(LinearModelOptions.af.opt)) {

            String option = aLine.getOptionValue(LinearModelOptions.af.opt);

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

        // The general covariates
        if (aLine.hasOption(LinearModelOptions.covariate_general.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.covariate_general.opt);

            covariatesGeneral = option.split(",");

        } else {

            covariatesGeneral = new String[0];

        }

        // The specific covariates
        if (aLine.hasOption(LinearModelOptions.covariate_specific.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.covariate_specific.opt);

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
        if (aLine.hasOption(LinearModelOptions.childId.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.childId.opt);

            if (option.length() == 0) {

                throw new IllegalArgumentException("Empty column name found for the child id.");

            }

            PhenotypesHandler.childIdColumn = option;

        }

        // The models
        if (aLine.hasOption(LinearModelOptions.model.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.model.opt);

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

            option = aLine.getOptionValue(LinearModelOptions.nVariants.opt);

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
        if (aLine.hasOption(LinearModelOptions.timeOut.opt)) {

            option = aLine.getOptionValue(LinearModelOptions.timeOut.opt);

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
        variantLog = aLine.hasOption(LinearModelOptions.variantLog.opt);

    }
}
