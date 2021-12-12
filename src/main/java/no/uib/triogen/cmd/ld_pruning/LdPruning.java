package no.uib.triogen.cmd.ld_pruning;

import java.io.IOException;
import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.processing.ld.pruning.SimpleLdPruner;

/**
 * Simple LD pruning of a result set.
 *
 * @author Marc Vaudel
 */
public class LdPruning {

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
            LdPruningOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            LdPruningOptionsBean bean = new LdPruningOptionsBean(commandLine);

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
     * @param bean the bean of command line parameters
     *
     * @throws IOException Exception thrown if an error occurs while reading the
     * ld file.
     */
    private static void run(
            LdPruningOptionsBean bean
    ) throws IOException {

        SimpleLdPruner pruner = new SimpleLdPruner(
                bean.ldMatrixFilePath,
                bean.resultsFile,
                bean.destinationFile,
                bean.minR2,
                bean.maxP,
                bean.pColNames,
                bean.idColName,
                bean.rsidColName,
                bean.phenoColName,
                bean.contigColName,
                bean.separator,
                bean.buildNumber,
                bean.ensemblPopulation,
                bean.ldLinkPopulation,
                bean.ldLinkToken
        );

        pruner.run();

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
            lPrintWriter.print("  Linkage Disequilibrium Pruning  " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The Linkage Disequilibrium Pruning command line prunes association results by LD." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(LdPruningOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
