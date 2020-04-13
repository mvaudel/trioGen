package no.uib.triogen.cmd.locus;

import no.uib.triogen.cmd.ld.*;
import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LocusZoomOptionsBean {

    /**
     * The ld matrix file.
     */
    public final File ldMatrixFile;
    /**
     * The results file.
     */
    public final File resultsFile;
    /**
     * The output file.
     */
    public final File outputFile;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public LocusZoomOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (LocusZoomOptions option : LocusZoomOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The ld matrix file
        String filePath = aLine.getOptionValue(LocusZoomOptions.ldMatrix.opt);

        ldMatrixFile = new File(filePath);

        if (!ldMatrixFile.exists()) {

            throw new IllegalArgumentException("Ld matrix file (" + ldMatrixFile + ") not found.");

        }

        // The ld matrix file
         filePath = aLine.getOptionValue(LocusZoomOptions.results.opt);

        resultsFile = new File(filePath);

        if (!resultsFile.exists()) {

            throw new IllegalArgumentException("Results file (" + resultsFile + ") not found.");

        }

        // The output file
         filePath = aLine.getOptionValue(LocusZoomOptions.out.opt);

        outputFile = new File(filePath);

        if (!outputFile.exists()) {

            throw new IllegalArgumentException("Output file (" + outputFile + ") not found.");

        }

    }
}
