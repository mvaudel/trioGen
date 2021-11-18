package no.uib.triogen.cmd.variant_file;

import java.io.File;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class VariantFileOptionsBean {

    /**
     * The path to the genotypes file with chromosome encoded.
     */
    public final String genotypesFilePath;
    /**
     * The file listing the variants to process.
     */
    public final File variantFile;
    /**
     * The build number.
     */
    public int buildNumber = 37;
    /**
     * The reference population to use for Ensembl.
     */
    public String ensemblPopulation = "1000GENOMES:phase_3:GBR";
    /**
     * The reference population to use for LDlink.
     */
    public String ldLinkPopulation = "CEU&GBR";
    /**
     * The token to use for LDlink.
     */
    public String ldLinkToken = null;
    /**
     * The source.
     */
    public final String source;
    /**
     * The minimal r2 for a proxy.
     */
    public double r2 = 0.2;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * The file where to write the missing ids.
     */
    public final File missingFile;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public VariantFileOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (VariantFileOptions option : VariantFileOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The genotype files
        genotypesFilePath = aLine.getOptionValue(VariantFileOptions.geno.opt);

        // The variant id file
        String filePath = aLine.getOptionValue(VariantFileOptions.variantId.opt);

        variantFile = new File(filePath);

        if (!variantFile.exists()) {

            throw new IllegalArgumentException("Variant id file (" + variantFile + ") not found.");

        }

        // The Ensembl build
        if (aLine.hasOption(VariantFileOptions.build.opt)) {

            String buildString = aLine.getOptionValue(VariantFileOptions.build.opt);

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
        if (aLine.hasOption(VariantFileOptions.ensemblPopulation.opt)) {

            ensemblPopulation = aLine.getOptionValue(VariantFileOptions.ensemblPopulation.opt);

        }

        // The reference population for proxies in LDlink
        if (aLine.hasOption(VariantFileOptions.ldlinkPopulation.opt)) {

            ldLinkPopulation = aLine.getOptionValue(VariantFileOptions.ldlinkPopulation.opt);

        }

        // The token for proxies in LDlink
        if (aLine.hasOption(VariantFileOptions.ldlinkToken.opt)) {

            ldLinkToken = aLine.getOptionValue(VariantFileOptions.ldlinkToken.opt);

        }

        // The minimal r2 for proxies
        if (aLine.hasOption(VariantFileOptions.r2.opt)) {

            String r2String = aLine.getOptionValue(VariantFileOptions.r2.opt);

            try {

                r2 = Double.parseDouble(r2String);

            } catch (Exception e) {

                throw new IllegalArgumentException("Input for r2 (" + r2String + ") cannot be parsed as a number.");

            }
        }

        // The output file
        filePath = aLine.getOptionValue(VariantFileOptions.out.opt);

        destinationFile = new File(filePath);

        File destinationFolder = destinationFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }

        missingFile = new File(destinationFolder, destinationFile.getName() + "_missing");

        // The source
        source = aLine.getOptionValue(VariantFileOptions.source.opt);

    }
}
