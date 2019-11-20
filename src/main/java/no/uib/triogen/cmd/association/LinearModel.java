package no.uib.triogen.cmd.association;

import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import static no.uib.triogen.io.Utils.lineSeparator;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.VariantList;
import no.uib.triogen.processing.association.linear_model.LinearModelComputer;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;

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

            run(bean);

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
            LinearModelOptionsBean bean
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);
        VariantList variantList = bean.variantFile == null ? null : VariantList.getVariantList(bean.variantFile);

        LinearModelComputer linearModelComputer = new LinearModelComputer(
                bean.genotypesFile,
                bean.genotypesFileType,
                variantList,
                childToParentMap,
                bean.phenotypesFile,
                bean.phenoNames,
                bean.destinationFile,
                bean.nVariants
        );

        try {

            linearModelComputer.run(
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

        try (PrintWriter lPrintWriter = new PrintWriter(System.out)) {
            lPrintWriter.print(lineSeparator);
            lPrintWriter.print("==================================" + lineSeparator);
            lPrintWriter.print("              trioGen             " + lineSeparator);
            lPrintWriter.print("               ****               " + lineSeparator);
            lPrintWriter.print("   Transmitted Allele Extraction  " + lineSeparator);
            lPrintWriter.print("==================================" + lineSeparator);
            lPrintWriter.print(lineSeparator
                    + "The extraction command line iterates a vcf file and extracts transmitted alleles of children." + lineSeparator
                    + lineSeparator
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + lineSeparator
                    + lineSeparator
                    + "----------------------"
                    + lineSeparator
                    + "OPTIONS"
                    + lineSeparator
                    + "----------------------" + lineSeparator
                    + lineSeparator);
            lPrintWriter.print(LinearModelOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
