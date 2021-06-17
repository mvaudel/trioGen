package no.uib.triogen.cmd.prs_train;

import java.io.File;
import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import no.uib.triogen.log.SimpleCliLogger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.processing.prs.PrsTrainer;

/**
 * Computes a simple polygenic risk score based on the sum of a list of betas.
 *
 * @author Marc Vaudel
 */
public class PrsTrain {

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
            PrsTrainOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            PrsTrainOptionsBean bean = new PrsTrainOptionsBean(commandLine);

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
     * @param command the command line as string
     */
    private static void run(
            PrsTrainOptionsBean bean,
            String command
    ) {

        String resultStem = bean.destinationFile.getAbsolutePath();

        if (resultStem.endsWith(".gz")) {

            resultStem = resultStem.substring(0, resultStem.length() - 3);

        }

        File logFile = new File(resultStem + ".log.gz");

        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "SimpleScore");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        try {

            PrsTrainer prsTrainer = new PrsTrainer(
                    bean.trainingFile,
                    bean.ldMatrixFilePath,
                    bean.destinationFile,
                    bean.snpIdColumn,
                    bean.chrColumn,
                    bean.posColumn,
                    bean.refColumn,
                    bean.eaColumn,
                    bean.betaPattern,
                    bean.sePattern,
                    bean.pPattern,
                    bean.model,
                    bean.variables,
                    bean.nSnpPerLocusThreshold,
                    bean.ldLocusThreshold,
                    bean.ldTopHitThreshold,
                    logger
            );

            prsTrainer.run();

        } catch (Throwable e) {

            e.printStackTrace();

        } finally {

            logger.close();

        }
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
            lPrintWriter.print("              PrsTrain            " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The PrsTrain command exports a list of weights from the pruning of trio summary statistics." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(PrsTrainOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
