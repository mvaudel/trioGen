package no.uib.triogen.cmd.ld_matrix;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import no.uib.triogen.TrioGen;
import no.uib.triogen.model.family.ChildToParentMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.processing.ld.LdMatrixComputer;

/**
 * Computes LD between variants and saves the results in a matrix.
 *
 * @author Marc Vaudel
 */
public class LdMatrix {

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
            LdMatrixOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            LdMatrixOptionsBean bean = new LdMatrixOptionsBean(commandLine);

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
     * @param bean The bean of command line parameters.
     * @param command The string version of the command.
     */
    public static void run(
            LdMatrixOptionsBean bean,
            String command
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(bean.chromosome);

        if (inheritanceMap == null) {

            throw new IllegalArgumentException("Mode of inheritance not implemented for " + bean.chromosome + ".");

        }
        
        int defaultMotherPlooidy = InheritanceUtils.getDefaultMotherPloidy(bean.chromosome);
        int defaultFatherPlooidy = InheritanceUtils.getDefaultFatherPloidy(bean.chromosome);

        File logFile = new File(bean.destinationFilePath + ".log.gz");
        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "LinearModel");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        LdMatrixComputer computer = new LdMatrixComputer(
                bean.genotypesFile,
                inheritanceMap,
                defaultMotherPlooidy,
                defaultFatherPlooidy,
                childToParentMap,
                bean.destinationFilePath,
                bean.maxDistance,
                bean.minR2,
                bean.nVariants,
                logger
        );

        try {

            computer.run(
                    bean.timeOut
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
            lPrintWriter.print("   Linkage Disequilibrium Matrix  " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The Linkage Disequilibrium Matrix command line iterates a vcf file and extracts variants in LD." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(LdMatrixOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
