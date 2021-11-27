package no.uib.triogen.cmd.locus;

import java.io.File;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class LocusZoomOptionsBean {

    /**
     * The phenotype of interest.
     */
    public String targetPhenotype = null;
    /**
     * The file listing the variants to process.
     */
    public File variantFile = null;
    /**
     * The maximum distance from the variant in bp.
     */
    public int maxDistance = 1000000;
    /**
     * The number of the build, e.g. 38 for GRCh38.
     */
    public int buildNumber = 38;
    /**
     * The ld matrix file.
     */
    public final File ldMatrixFile;
    /**
     * The results file.
     */
    public final File resultsFile;
    /**
     * Stem for the output file.
     */
    public final String outputFileStem;
    /**
     * Stem for the gene coordinates file.
     */
    public String geneCoordinatesFileStem = null;
    /**
     * Path for the file where to write the log.
     */
    public String logFilePath = null;

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

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The phenotype
        if (CliUtils.hasOption(aLine, LocusZoomOptions.phenoName)) {

            targetPhenotype = CliUtils.getOptionValue(aLine, LocusZoomOptions.phenoName);

        }

        // The variant ids
        if (CliUtils.hasOption(aLine, LocusZoomOptions.variantId)) {

            String filePath = CliUtils.getOptionValue(aLine, LocusZoomOptions.variantId);

            variantFile = new File(filePath);

            if (!variantFile.exists()) {

                throw new IllegalArgumentException("Variant file (" + variantFile + ") not found.");

            }
        }

        // The max distance
        if (CliUtils.hasOption(aLine, LocusZoomOptions.maxDistance)) {

            String stringValue = CliUtils.getOptionValue(aLine, LocusZoomOptions.maxDistance);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance <= 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a stricly positive integer.");

            }
        }

        // The build number
        if (CliUtils.hasOption(aLine, LocusZoomOptions.buildNumber)) {

            String stringValue = CliUtils.getOptionValue(aLine, LocusZoomOptions.buildNumber);

            buildNumber = Integer.parseInt(stringValue);

            if (buildNumber != 37 && buildNumber != 38) {

                throw new IllegalArgumentException("Build (" + buildNumber + ") not supported, must be 37 or 38.");

            }
        }

        // The ld matrix file
        String filePath = CliUtils.getOptionValue(aLine, LocusZoomOptions.ldMatrix);

        ldMatrixFile = new File(filePath);

        if (!ldMatrixFile.exists()) {

            throw new IllegalArgumentException("Ld matrix file (" + ldMatrixFile + ") not found.");

        }

        // The association results file
        filePath = CliUtils.getOptionValue(aLine, LocusZoomOptions.results);

        resultsFile = new File(filePath);

        if (!resultsFile.exists()) {

            throw new IllegalArgumentException("Results file (" + resultsFile + ") not found.");

        }

        // The output file
        outputFileStem = CliUtils.getOptionValue(aLine, LocusZoomOptions.out);

        // The gene coordinates file
        if (CliUtils.hasOption(aLine, LocusZoomOptions.geneCoordinates)) {

            geneCoordinatesFileStem = CliUtils.getOptionValue(aLine, LocusZoomOptions.geneCoordinates);
            
        }

        // The log file
        if (CliUtils.hasOption(aLine, LocusZoomOptions.log)) {

            logFilePath = CliUtils.getOptionValue(aLine, LocusZoomOptions.log);
            
        }
    }
}
