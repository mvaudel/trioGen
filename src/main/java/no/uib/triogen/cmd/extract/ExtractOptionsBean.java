package no.uib.triogen.cmd.extract;

import java.io.File;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class ExtractOptionsBean {

    /**
     * The results file to process.
     */
    public final File inputFile;
    /**
     * Splits the results by variant.
     */
    public final boolean splitByVariant;
    /**
     * Splits the results by phenotype.
     */
    public final boolean splitByPheno;
    /**
     * The columns.
     */
    public final String[] columns;
    /**
     * The variant ids.
     */
    public final String[] variantIds;
    /**
     * The phenotype names.
     */
    public final String[] phenoNames;
    /**
     * The stem of the output file.
     */
    public final String outputStem;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public ExtractOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (ExtractOptions option : ExtractOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The input file
        String filePath = CliUtils.getOptionValue(aLine, ExtractOptions.input);

        inputFile = new File(filePath);

        if (!inputFile.exists()) {

            throw new IllegalArgumentException("Input file (" + inputFile + ") not found.");

        }

        // The output stem
        filePath = CliUtils.getOptionValue(aLine, ExtractOptions.output);

        outputStem = filePath;

        File destinationFolder = (new File(outputStem)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Split by variant
        splitByVariant = CliUtils.hasOption(aLine, ExtractOptions.split_by_variant);

        // Split by pheno
        splitByPheno = CliUtils.hasOption(aLine, ExtractOptions.split_by_pheno);

        // The values columns
        if (CliUtils.hasOption(aLine, ExtractOptions.columns)) {

            String option = CliUtils.getOptionValue(aLine, ExtractOptions.columns);
            
            columns = option.split(",");
            
        } else {
            
            columns = null;
            
        }

        // The variant ids
        if (CliUtils.hasOption(aLine, ExtractOptions.variantId)) {

            String option = CliUtils.getOptionValue(aLine, ExtractOptions.variantId);
            
            variantIds = option.split(",");
            
        } else {
            
            variantIds = null;
            
        }

        // The pheno names
        if (CliUtils.hasOption(aLine, ExtractOptions.phenoName)) {

            String option = CliUtils.getOptionValue(aLine, ExtractOptions.phenoName);
            
            phenoNames = option.split(",");
            
        } else {
            
            phenoNames = null;
            
        }
    }
}
