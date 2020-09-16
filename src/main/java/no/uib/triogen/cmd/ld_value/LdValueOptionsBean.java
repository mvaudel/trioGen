package no.uib.triogen.cmd.ld_value;

import java.io.File;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LdValueOptionsBean {

    /**
     * The ld matrix file path.
     */
    public final String ldMatrixFilePath;
    /**
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * File where to write the output.
     */
    public final File destinationFile;
    /**
     * The minimal r2 to export.
     */
    public Double minR2 = null;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public LdValueOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (LdValueOptions option : LdValueOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The ld matrix file
        String filePath = aLine.getOptionValue(LdValueOptions.ldMatrix.opt);

        ldMatrixFilePath = filePath;

        // The variant ids
        if (aLine.hasOption(LdValueOptions.variantId.opt)) {

            filePath = aLine.getOptionValue(LdValueOptions.variantId.opt);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The output file
        filePath = aLine.getOptionValue(LdValueOptions.out.opt);
        
        if (!filePath.endsWith(".gz")) {
            
            filePath = filePath + ".gz";
            
        }

        destinationFile = new File(filePath);

        File destinationFolder = destinationFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }
        
        // The min r2 to report
        if (aLine.hasOption(LdValueOptions.minR2.opt)) {
            
            String stringValue = aLine.getOptionValue(LdValueOptions.minR2.opt);
            
            try {
                
                minR2 = Double.valueOf(stringValue);
                
            } catch (Exception e) {
                
            throw new IllegalArgumentException("Input for minimal r2 cannot be parsed as a number (" + stringValue + ").");
                
            }
        }
    }
}
