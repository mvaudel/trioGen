package no.uib.triogen.cmd.results;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.model.geno.Model;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import no.uib.triogen.processing.association.linear_model.LinearModelRunnable;
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
     * The category columns.
     */
    public final String[] categoryColumns;
    /**
     * The value columns.
     */
    public final String[] valueColumns;
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

        // The category columns
        if (aLine.hasOption(ExtractOptions.category.opt)) {

            String option = aLine.getOptionValue(ExtractOptions.category.opt);
            
            categoryColumns = option.split(",");
            
        } else {
            
            categoryColumns = null;
            
        }

        // The values columns
        if (aLine.hasOption(ExtractOptions.value.opt)) {

            String option = aLine.getOptionValue(ExtractOptions.value.opt);
            
            valueColumns = option.split(",");
            
        } else {
            
            valueColumns = null;
            
        }

    }
}
