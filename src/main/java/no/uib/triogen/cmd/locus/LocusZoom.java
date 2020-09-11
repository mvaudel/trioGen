package no.uib.triogen.cmd.locus;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import no.uib.triogen.export.LocusZoomExtractor;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * Extracts association data and ld for a locus of interest.
 *
 * @author Marc Vaudel
 */
public class LocusZoom {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        if (args.length == 0
                || args.length == 1 && args[0].equals("-h")
                || args.length == 1 && args[0].equals("--help")) {

            printHelp();
            return;

        }

        if (args.length == 1 && args[0].equals("-v")
                || args.length == 1 && args[0].equals("--version")) {

            System.out.println(TrioGen.getVersion());

            return;

        }

        try {

            Options lOptions = new Options();
            LocusZoomOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            LocusZoomOptionsBean bean = new LocusZoomOptionsBean(commandLine);

            run(
                    bean
            );

        } catch (Throwable e) {

            e.printStackTrace();
        }
    }

    /**
     * Runs the command.
     *
     * @param bean The bean of command line parameters.
     *
     * @throws IOException Exception thrown if an error occurred while reading
     * or writing a file.
     */
    private static void run(
            LocusZoomOptionsBean bean
    ) throws IOException {

        File logFile = bean.logFilePath == null ? new File(bean.outputFileStem + ".log") : new File(bean.logFilePath);

        SimpleCliLogger logger = new SimpleCliLogger(logFile);

        VariantList variantList = VariantList.getVariantList(bean.variantFile);

        try {

            LocusZoomExtractor.writeData(
                    bean.targetPhenotype,
                    variantList,
                    bean.maxDistance,
                    bean.buildNumber,
                    bean.resultsFile,
                    bean.ldMatrixFile,
                    bean.outputFileStem,
                    bean.geneCoordinatesFileStem,
                    logger
            );

        } catch (Throwable t) {

            t.printStackTrace();

            logger.logError(t.getLocalizedMessage());

        }

        logger.close();

    }

    /**
     * Prints basic help
     */
    private static void printHelp() {

        try ( PrintWriter lPrintWriter = new PrintWriter(System.out)) {
            lPrintWriter.print(LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print("              trioGen             " + LINE_SEPARATOR);
            lPrintWriter.print("               ****               " + LINE_SEPARATOR);
            lPrintWriter.print("             Locus Zoom           " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The Locus Zoom command line extracts association data and ld for a locus of interest." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(LocusZoomOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
