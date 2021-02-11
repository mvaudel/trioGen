package no.uib.triogen.cmd.vcf_to_bgen;

import java.io.File;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class VcfToBgenOptionsBean {

    /**
     * The vcf file.
     */
    public final File vcfFile;
    /**
     * The bgen file.
     */
    public final File bgenFile;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public VcfToBgenOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (VcfToBgenOptions option : VcfToBgenOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The input file
        String filePath = aLine.getOptionValue(VcfToBgenOptions.input.opt);

        vcfFile = new File(filePath);

        if (!vcfFile.exists()) {

            throw new IllegalArgumentException("Vcf file (" + vcfFile + ") not found.");

        }

        // The output stem
        filePath = aLine.getOptionValue(VcfToBgenOptions.output.opt);

        bgenFile = new File(filePath);

        File destinationFolder = bgenFile.getParentFile();

        if (!destinationFolder.exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFolder + ") not found.");

        }
    }
}
