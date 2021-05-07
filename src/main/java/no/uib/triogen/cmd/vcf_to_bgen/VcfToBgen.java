package no.uib.triogen.cmd.vcf_to_bgen;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import no.uib.triogen.TrioGen;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.conversion.VcfToBgenConverter;

/**
 * Converts a vcf file to a bgen file.
 *
 * @author Marc Vaudel
 */
public class VcfToBgen {

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
            VcfToBgenOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            VcfToBgenOptionsBean bean = new VcfToBgenOptionsBean(commandLine);

            run(bean);

        } catch (Throwable e) {

            e.printStackTrace();
        }
    }

    /**
     * Runs the command.
     *
     * @param bean the bean of command line parameters
     *
     * @throws IOException Exception thrown if an I/O error occurs.
     */
    private static void run(
            VcfToBgenOptionsBean bean
    ) throws IOException {

        VcfToBgenConverter vcfToBgenConverter = new VcfToBgenConverter();
        vcfToBgenConverter.convert(bean.vcfFile, bean.bgenFile);

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
            lPrintWriter.print("      Vcf to Bgen Conversion      " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The VcfToBgen command line writes the phased genotypes from a vcf to the bgen v.1.3 format with zstd compression and one probability per allele." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(VcfToBgenOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }

    /**
     * Returns the default columns that are always included in the results.
     *
     * @return The default columns that are always included in the results.
     */
    public static HashSet<String> getDefualtColumns() {

        HashSet<String> results = new HashSet<>(4);
        results.add("phenotype");
        results.add("variantId");
        results.add("n");
        results.add("nAlt");
        results.add("nH");

        return results;

    }
}
