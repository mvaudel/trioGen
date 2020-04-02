package no.uib.triogen.cmd.association;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.log.Logger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.Model;
import no.uib.triogen.model.geno.VariantList;
import no.uib.triogen.processing.association.linear_model.LinearModelComputer;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;

/**
 * Runs multiple linear models for the association with phenotypes.
 *
 * @author Marc Vaudel
 */
public class LinearModel {

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
            LinearModelOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            LinearModelOptionsBean bean = new LinearModelOptionsBean(commandLine);

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
            LinearModelOptionsBean bean,
            String command
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);
        VariantList variantList = bean.variantFile == null ? null : VariantList.getVariantList(bean.variantFile);
        Model[] models = Arrays.stream(bean.modelNames)
                .map(
                        modelName -> Model.valueOf(modelName)
                )
                .toArray(Model[]::new);

        String resultStem = bean.destinationFile.getAbsolutePath();

        if (resultStem.endsWith(".gz")) {

            resultStem = resultStem.substring(0, resultStem.length() - 3);

        }

        File logFile = new File(resultStem + ".log.gz");
        File variantLogFile = bean.variantLog ? new File(resultStem + ".variantLog.gz") : null;

        Logger logger = new Logger(logFile, variantLogFile);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "LinearModel");
        logger.writeComment("Arguments", command);
        logger.writeHeaders();

        LinearModelComputer linearModelComputer = new LinearModelComputer(
                bean.genotypesFile,
                bean.genotypesFileType,
                variantList,
                bean.maf,
                childToParentMap,
                bean.phenotypesFile,
                bean.phenoNames,
                bean.covariates,
                models,
                bean.destinationFile,
                bean.nVariants,
                logger
        );

        try {

            linearModelComputer.run(
                    bean.timeOut,
                    bean.test
            );

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
            lPrintWriter.print("      Linear Model Regression     " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The linear model regression command performs linear regression using various models." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(LinearModelOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
