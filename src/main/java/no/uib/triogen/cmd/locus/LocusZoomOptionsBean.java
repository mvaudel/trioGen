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
     * The phenotype of interest.
     */
    public final String targetPhenotype;
    /**
     * The id of the variant of interest.
     */
    public final String targetVariant;
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
     * The output file.
     */
    public final File outputFile;
    /**
     * The gene coordinates file.
     */
    public File geneCoordinatesFile = null;

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

        // The phenotype
        targetPhenotype = aLine.getOptionValue(LocusZoomOptions.phenoName.opt);

        // The variant id
        targetVariant = aLine.getOptionValue(LocusZoomOptions.variantId.opt);

        // The max distance
        if (aLine.hasOption(LocusZoomOptions.maxDistance.opt)) {

            String stringValue = aLine.getOptionValue(LocusZoomOptions.maxDistance.opt);

            maxDistance = Integer.parseInt(stringValue);

            if (maxDistance <= 0) {

                throw new IllegalArgumentException("Distance (" + maxDistance + ") must be a stricly positive integer.");

            }
        }

        // The build number
        if (aLine.hasOption(LocusZoomOptions.buildNumber.opt)) {

            String stringValue = aLine.getOptionValue(LocusZoomOptions.buildNumber.opt);

            buildNumber = Integer.parseInt(stringValue);

            if (buildNumber != 37 && buildNumber != 38) {

                throw new IllegalArgumentException("Build (" + buildNumber + ") not supported, must be 37 or 38.");

            }
        }

        // The ld matrix file
        String filePath = aLine.getOptionValue(LocusZoomOptions.ldMatrix.opt);

        ldMatrixFile = new File(filePath);

        if (!ldMatrixFile.exists()) {

            throw new IllegalArgumentException("Ld matrix file (" + ldMatrixFile + ") not found.");

        }

        // The association results file
        filePath = aLine.getOptionValue(LocusZoomOptions.results.opt);

        resultsFile = new File(filePath);

        if (!resultsFile.exists()) {

            throw new IllegalArgumentException("Results file (" + resultsFile + ") not found.");

        }

        // The output file
        filePath = aLine.getOptionValue(LocusZoomOptions.out.opt);

        outputFile = new File(filePath);

        if (!outputFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + outputFile.getParentFile() + ") not found.");

        }

        // The gene coordinates file
        if (aLine.hasOption(LocusZoomOptions.geneCoordinates.opt)) {

            filePath = aLine.getOptionValue(LocusZoomOptions.geneCoordinates.opt);

            geneCoordinatesFile = new File(filePath);

            if (!geneCoordinatesFile.getParentFile().exists()) {

                throw new IllegalArgumentException("Gene coordinates folder (" + geneCoordinatesFile.getParentFile() + ") not found.");

            }
        }
    }
}
