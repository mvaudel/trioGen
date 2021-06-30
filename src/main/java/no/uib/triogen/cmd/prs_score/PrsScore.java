package no.uib.triogen.cmd.prs_score;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import no.uib.triogen.TrioGen;
import no.uib.triogen.log.SimpleCliLogger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.processing.prs.PrsScorer;
import no.uib.triogen.processing.prs.PrsTrainer;

/**
 * Computes a simple polygenic risk score based on the sum of a list of betas.
 *
 * @author Marc Vaudel
 */
public class PrsScore {

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
            PrsScoreOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            PrsScoreOptionsBean bean = new PrsScoreOptionsBean(commandLine);

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
            PrsScoreOptionsBean bean,
            String command
    ) {

        String filePath = bean.destinationFile.getAbsolutePath();
        String logPath;

        if (!filePath.endsWith(".gz")) {

            logPath = filePath;
            filePath = filePath + ".gz";

        } else {
            
            logPath = filePath.substring(0, filePath.length() - 3);
            
        }

        File destinationFile = new File(filePath);
        File logFile = new File(logPath + ".log.gz");

        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "PrsScore");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        VariantList variantList = null;

        if (bean.variantFile != null) {

            variantList = VariantList.getVariantList(
                    bean.variantFile,
                    null
            );

            if (variantList.variantId.length == 0) {

                logger.logMessage("No target variant in file " + bean.variantFile + ".");

            }

            variantList.index(0);

        }

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        try {

            PrsScorer prsScorer = new PrsScorer(
                    bean.genotypesFile,
                    childToParentMap,
                    bean.scoreFile,
                    destinationFile,
                    variantList,
                    bean.betaPattern,
                    bean.sePattern,
                    bean.model,
                    bean.variables,
                    bean.pValueThreshold,
                    bean.betaQuantile,
                    bean.scoringMode,
                    logger
            );

            prsScorer.run();

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
            lPrintWriter.print("              PrsScore            " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The PrsScore command computes a score on genotypes based on a PrsTrain file." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(PrsScoreOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
