package no.uib.triogen.cmd.locus;

import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;

/**
 * Computes LD between variants and saves the results in a matrix.
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
                    bean,
                    String.join(" ", args)
                    );

        } catch (Throwable e) {

            e.printStackTrace();
        }
    }

    /**
     * Runs the command.
     *
     * @param bean the bean of command line parameters
     */
    private static void run(
            LocusZoomOptionsBean bean,
            String command
    ) {
//
//        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);
//
//        File logFile = new File(bean.destinationFilePath + ".log.gz");
//        Logger logger = new Logger(logFile, null);
//        logger.writeComment("Software", "TrioGen");
//        logger.writeComment("Version", TrioGen.getVersion());
//        logger.writeComment("Command", "LinearModel");
//        logger.writeComment("Arguments", command);
//        logger.writeHeaders();
//        
//        File outputFile = new File(bean.destinationFilePath + ".tld");
//
//        LdMatrixComputer computer = new LdMatrixComputer(
//                bean.genotypesFile, 
//                bean.genotypesFileType, 
//                childToParentMap, 
//                outputFile, 
//                bean.maxDistance, 
//                bean.minR2,
//                bean.maf,
//                bean.hardCalls, 
//                bean.nVariants, 
//                bean.downstreamLoadingFactor,
//                bean.upstreamLoadingFactor,
//                logger
//        );
//
//        try {
//
//            computer.run(
//                    bean.timeOut,
//                    bean.testIteration,
//                    bean.test
//            );
//
//        } catch (Throwable e) {
//
//            e.printStackTrace();
//
//        }
    }

    /**
     * Prints basic help
     */
    private static void printHelp() {

        try (PrintWriter lPrintWriter = new PrintWriter(System.out)) {
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
