package no.uib.triogen.cmd.transmission;

import java.io.PrintWriter;
import no.uib.triogen.TrioGen;
import static no.uib.triogen.io.Utils.lineSeparator;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.transmission.extraction.Extractor;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;

/**
 * Extracts h matrices from the vcf files.
 *
 * @author Marc Vaudel
 */
public class ExtractTransmission {

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
            ExtractTransmissionOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            ExtractTransmissionOptionsBean bean = new ExtractTransmissionOptionsBean(commandLine);
            
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
            ExtractTransmissionOptionsBean bean
    ) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        Extractor extractor = new Extractor(
                bean.genotypesFile, 
                bean.genotypesFileType,
                childToParentMap, 
                bean.destinationStem,
                bean.nVariants
        );

        try {

            extractor.run(
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
            lPrintWriter.print(ExtractTransmissionOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
