package no.uib.triogen.cmd.mendelian_check;

import java.io.File;
import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.processing.mendelian_check.MendelianCheckComputer;

/**
 * Writes a report on Mendelian errors in the given data set.
 *
 * @author Marc Vaudel
 */
public class MendelianCheck {

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
            MendelianCheckOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            MendelianCheckOptionsBean bean = new MendelianCheckOptionsBean(commandLine);

            if (bean.test) {

                CustomVcfIterator.nLimit = 1000;

            }

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
            MendelianCheckOptionsBean bean,
            String command
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        File logFile = new File(bean.destinationFile.getAbsolutePath() + ".log.gz");
        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "LinearModel");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();
        
        VariantList variantList = bean.variantFile == null ? null : VariantList.getVariantList(bean.variantFile);

        MendelianCheckComputer computer = new MendelianCheckComputer(
                bean.genotypesFile,
                bean.genotypesFileType,
                variantList,
                childToParentMap,
                bean.destinationFile,
                bean.maf,
                bean.nVariants,
                logger
        );

        try {

            computer.run(
                    bean.timeOut,
                    bean.test
            );

        } catch (Throwable e) {

            e.printStackTrace();

        }
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
            lPrintWriter.print("          Mendelian Check         " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The Mendelian check command line iterates a vcf file and extracts summary statistics on Mendelian errors." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(MendelianCheckOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
