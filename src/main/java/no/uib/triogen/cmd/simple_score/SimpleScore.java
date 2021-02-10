package no.uib.triogen.cmd.simple_score;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.Collectors;
import no.uib.triogen.TrioGen;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.model.simple_score.VariantWeightList;
import no.uib.triogen.processing.simple_score.SimpleScoreComputer;

/**
 * Runs multiple linear models for the association with phenotypes.
 *
 * @author Marc Vaudel
 */
public class SimpleScore {

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
            SimpleScoreOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            SimpleScoreOptionsBean bean = new SimpleScoreOptionsBean(commandLine);

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
            SimpleScoreOptionsBean bean,
            String command
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);
        VariantWeightList variantWeightList = VariantWeightList.getVariantWeightList(bean.variantFile);

        HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(bean.chromosome);

        if (inheritanceMap == null) {

            throw new IllegalArgumentException("Mode of inheritance not implemented for " + bean.chromosome + ".");

        }
        
        int defaultMotherPlooidy = InheritanceUtils.getDefaultMotherPloidy(bean.chromosome);
        int defaultFatherPlooidy = InheritanceUtils.getDefaultFatherPloidy(bean.chromosome);

        String resultStem = bean.destinationFile.getAbsolutePath();

        if (resultStem.endsWith(".gz")) {

            resultStem = resultStem.substring(0, resultStem.length() - 3);

        }

        File logFile = new File(resultStem + ".log.gz");
        File variantLogFile = bean.variantLog ? new File(resultStem + ".variantLog.gz") : null;

        SimpleCliLogger logger = new SimpleCliLogger(logFile, variantLogFile);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "SimpleScore");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        SimpleScoreComputer scoreComputer = new SimpleScoreComputer(
                bean.genotypesFile, 
                inheritanceMap, 
                defaultMotherPlooidy,
                defaultFatherPlooidy,
                childToParentMap, 
                variantWeightList, 
                bean.phenotypesFile, 
                bean.phenoNames, 
                bean.destinationFile, 
                logger
        );

        try {

            scoreComputer.computeScore();

        } catch (Throwable e) {

            logger.logError(
                    Arrays.stream(e.getStackTrace())
                            .map(
                                    element -> element.toString()
                            )
                            .collect(Collectors.joining(" "))
            );

            e.printStackTrace();

        }

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
            lPrintWriter.print("            Simple Score          " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The simple score command computes a simple risk score in trios based on a list of weights." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(SimpleScoreOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
