package no.uib.triogen.cmd.variant_file;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import no.uib.triogen.log.SimpleCliLogger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.variant_id.VariantListToFile;

/**
 * Makes a variant file from a list of rsids.
 *
 * @author Marc Vaudel
 */
public class VariantFile {

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
            VariantFileOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            VariantFileOptionsBean bean = new VariantFileOptionsBean(commandLine);

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
     * 
     * @throws IOException Exception thrown if an error occurred while reading or writing a file.
     */
    private static void run(
            VariantFileOptionsBean bean,
            String command
    ) throws IOException {

        File logFile = new File(bean.destinationFile.getParentFile(), bean.destinationFile.getName() + ".log.gz");

        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "VariantFile");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        VariantListToFile variantListToFile = new VariantListToFile();
        
        variantListToFile.writeVariantFile(
                bean.source, 
                bean.genotypesFilePath, 
                bean.variantFile, 
                bean.destinationFile, 
                bean.missingFile, 
                bean.buildNumber, 
                bean.ldLinkPopulation,
                bean.ldLinkToken,
                bean.ensemblPopulation, 
                bean.r2, 
                logger
        );

        logger.close();

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
            lPrintWriter.print("            Variant File          " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The variant file command creates a target file that can be used in other command lines from a list of variants." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(VariantFileOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
