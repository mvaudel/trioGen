package no.uib.triogen.cmd.results;

import java.io.File;
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

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The input file
        String filePath = aLine.getOptionValue(ExtractOptions.input.opt);

        inputFile = new File(filePath);

        if (!inputFile.exists()) {

            throw new IllegalArgumentException("Input file (" + inputFile + ") not found.");

        }

        // The output stem
        filePath = aLine.getOptionValue(ExtractOptions.output.opt);

        outputStem = filePath;

        File destinationFolder = (new File(outputStem)).getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        // Split by variant
        splitByVariant = aLine.hasOption(ExtractOptions.split_by_variant.opt);

        // Split by pheno
        splitByPheno = aLine.hasOption(ExtractOptions.split_by_pheno.opt);

        // The values columns
        if (aLine.hasOption(ExtractOptions.columns.opt)) {

            String option = aLine.getOptionValue(ExtractOptions.columns.opt);
            
            columns = option.split(",");
            
        } else {
            
            columns = null;
            
        }

        // The variant ids
        if (aLine.hasOption(ExtractOptions.variantId.opt)) {

            String option = aLine.getOptionValue(ExtractOptions.variantId.opt);
            
            variantIds = option.split(",");
            
        } else {
            
            variantIds = null;
            
        }

        // The pheno names
        if (aLine.hasOption(ExtractOptions.phenoName.opt)) {

            String option = aLine.getOptionValue(ExtractOptions.phenoName.opt);
            
            phenoNames = option.split(",");
            
        } else {
            
            phenoNames = null;
            
        }
    }
}
