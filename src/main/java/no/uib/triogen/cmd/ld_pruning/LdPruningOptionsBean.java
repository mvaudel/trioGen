package no.uib.triogen.cmd.ld_pruning;

import java.io.File;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LdPruningOptionsBean {

    /**
     * The ld matrix file path.
     */
    public final String ldMatrixFilePath;
    /**
     * The results file to prune.
     */
    public final File resultsFile;
    /**
     * File where to write the output.
     */
    public final File destinationFile;
    /**
     * The minimal r2 to export.
     */
    public double minR2 = 0.05;
    /**
     * The maximal p to export.
     */
    public double maxP = 1e-6;
    /**
     * The name of the variant id column.
     */
    public String idColName = "variantId";
    /**
     * The name of the p-value column.
     */
    public String pColName = "h.intercept.p";
    /**
     * The name of the pheno column.
     */
    public String phenoColName = null;
    /**
     * The name of the contig column.
     */
    public String contigColName = "contig";
    /**
     * The column separator.
     */
    public String separator = "\t";
    
    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public LdPruningOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (LdPruningOptions option : LdPruningOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The results file
        String filePath = aLine.getOptionValue(LdPruningOptions.results.opt);

        resultsFile = new File(filePath);
        
        if (!resultsFile.exists()) {
            
            throw new IllegalArgumentException("Results file " + resultsFile.getAbsolutePath() + " not found.");
            
        }
        

        // The ld matrix file
        filePath = aLine.getOptionValue(LdPruningOptions.ldMatrix.opt);

        ldMatrixFilePath = filePath;
        

        // The output file
        filePath = aLine.getOptionValue(LdPruningOptions.out.opt);
        
        if (!filePath.endsWith(".gz")) {
            
            filePath = filePath + ".gz";
            
        }

        destinationFile = new File(filePath);

        File destinationFolder = destinationFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }
        
        
        // The minimal r2 to consider in ld
        if (aLine.hasOption(LdPruningOptions.minR2.opt)) {
            
            String stringValue = aLine.getOptionValue(LdPruningOptions.minR2.opt);
            
            try {
                
                minR2 = Double.valueOf(stringValue);
                
            } catch (Exception e) {
                
            throw new IllegalArgumentException("Input for minimal r2 cannot be parsed as a number (" + stringValue + ").");
                
            }
        }
        
        
        // The maximal p to report
        if (aLine.hasOption(LdPruningOptions.maxP.opt)) {
            
            String stringValue = aLine.getOptionValue(LdPruningOptions.maxP.opt);
            
            try {
                
                maxP = Double.valueOf(stringValue);
                
            } catch (Exception e) {
                
            throw new IllegalArgumentException("Input for maximal p cannot be parsed as a number (" + stringValue + ").");
                
            }
        }
        
        
        // The name of the variant id column
        if (aLine.hasOption(LdPruningOptions.idColName.opt)) {
            
            idColName = aLine.getOptionValue(LdPruningOptions.idColName.opt);
            
        }
        
        
        // The name of the p-value column
        if (aLine.hasOption(LdPruningOptions.pColName.opt)) {
            
            pColName = aLine.getOptionValue(LdPruningOptions.pColName.opt);
            
        }
        
        
        // The name of the pheno column
        if (aLine.hasOption(LdPruningOptions.phenoColName.opt)) {
            
            phenoColName = aLine.getOptionValue(LdPruningOptions.phenoColName.opt);
            
        }
        
        
        // The name of the pheno column
        if (aLine.hasOption(LdPruningOptions.contigColName.opt)) {
            
            contigColName = aLine.getOptionValue(LdPruningOptions.contigColName.opt);
            
        }
        
        
        // The name of the pheno column
        if (aLine.hasOption(LdPruningOptions.separator.opt)) {
            
            separator = aLine.getOptionValue(LdPruningOptions.separator.opt);
            
        }
    }
}
