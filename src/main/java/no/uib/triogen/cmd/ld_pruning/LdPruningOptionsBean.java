package no.uib.triogen.cmd.ld_pruning;

import java.io.File;
import no.uib.triogen.utils.cli.CliUtils;
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
    public String ldMatrixFilePath = null;
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
     * The name of the variant rsid column.
     */
    public String rsidColName = "variantId";
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
     * The build number.
     */
    public int buildNumber = 37;
    /**
     * The reference population to use for Ensembl.
     */
    public String ensemblPopulation = null;
    /**
     * The reference population to use for LDlink.
     */
    public String ldLinkPopulation = null;
    /**
     * The token to use for LDlink.
     */
    public String ldLinkToken = null;
    
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

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The results file
        String filePath = CliUtils.getOptionValue(aLine, LdPruningOptions.results);

        resultsFile = new File(filePath);
        
        if (!resultsFile.exists()) {
            
            throw new IllegalArgumentException("Results file " + resultsFile.getAbsolutePath() + " not found.");
            
        }
        

        // The ld matrix file
        filePath = CliUtils.getOptionValue(aLine, LdPruningOptions.ldMatrix);

        ldMatrixFilePath = filePath;
        

        // The output file
        filePath = CliUtils.getOptionValue(aLine, LdPruningOptions.out);
        
        if (!filePath.endsWith(".gz")) {
            
            filePath = filePath + ".gz";
            
        }

        destinationFile = new File(filePath);

        File destinationFolder = destinationFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }
        
        
        // The minimal r2 to consider in ld
        if (CliUtils.hasOption(aLine, LdPruningOptions.minR2)) {
            
            String stringValue = CliUtils.getOptionValue(aLine, LdPruningOptions.minR2);
            
            try {
                
                minR2 = Double.valueOf(stringValue);
                
            } catch (Exception e) {
                
            throw new IllegalArgumentException("Input for minimal r2 cannot be parsed as a number (" + stringValue + ").");
                
            }
        }
        
        
        // The maximal p to report
        if (CliUtils.hasOption(aLine, LdPruningOptions.maxP)) {
            
            String stringValue = CliUtils.getOptionValue(aLine, LdPruningOptions.maxP);
            
            try {
                
                maxP = Double.valueOf(stringValue);
                
            } catch (Exception e) {
                
            throw new IllegalArgumentException("Input for maximal p cannot be parsed as a number (" + stringValue + ").");
                
            }
        }
        
        
        // The name of the variant id column
        if (CliUtils.hasOption(aLine, LdPruningOptions.idColName)) {
            
            idColName = CliUtils.getOptionValue(aLine, LdPruningOptions.idColName);
            
        }
        
        
        // The name of the rsid column
        if (CliUtils.hasOption(aLine, LdPruningOptions.rsidColName)) {
            
            rsidColName = CliUtils.getOptionValue(aLine, LdPruningOptions.rsidColName);
            
        }
        
        
        // The name of the p-value column
        if (CliUtils.hasOption(aLine, LdPruningOptions.pColName)) {
            
            pColName = CliUtils.getOptionValue(aLine, LdPruningOptions.pColName);
            
        }
        
        
        // The name of the pheno column
        if (CliUtils.hasOption(aLine, LdPruningOptions.phenoColName)) {
            
            phenoColName = CliUtils.getOptionValue(aLine, LdPruningOptions.phenoColName);
            
        }
        
        
        // The name of the pheno column
        if (CliUtils.hasOption(aLine, LdPruningOptions.contigColName)) {
            
            contigColName = CliUtils.getOptionValue(aLine, LdPruningOptions.contigColName);
            
        }
        
        
        // The name of the pheno column
        if (CliUtils.hasOption(aLine, LdPruningOptions.separator)) {
            
            separator = CliUtils.getOptionValue(aLine, LdPruningOptions.separator);
            
        }

        // The Ensembl build
        if (CliUtils.hasOption(aLine, LdPruningOptions.build)) {

            String buildString = CliUtils.getOptionValue(aLine, LdPruningOptions.build);

            try {

                buildNumber = Integer.parseInt(buildString);
                
                if (buildNumber != 37 && buildNumber != 38) {
                    
                throw new IllegalArgumentException("Input for build number (" + buildNumber + ") not supported, only 37 and 38 are currently supported.");
                    
                }

            } catch (Exception e) {

                throw new IllegalArgumentException("Input for build number (" + buildString + ") cannot be parsed as a number.");

            }
        }

        // The reference population for proxies in Ensembl
        if (CliUtils.hasOption(aLine, LdPruningOptions.ensemblPopulation)) {

            ensemblPopulation = CliUtils.getOptionValue(aLine, LdPruningOptions.ensemblPopulation);

        }

        // The reference population for proxies in LDlink
        if (CliUtils.hasOption(aLine, LdPruningOptions.ldlinkPopulation)) {

            ldLinkPopulation = CliUtils.getOptionValue(aLine, LdPruningOptions.ldlinkPopulation);

        }

        // The token for proxies in LDlink
        if (CliUtils.hasOption(aLine, LdPruningOptions.ldlinkToken)) {

            ldLinkToken = CliUtils.getOptionValue(aLine, LdPruningOptions.ldlinkToken);

        }
    }
}
