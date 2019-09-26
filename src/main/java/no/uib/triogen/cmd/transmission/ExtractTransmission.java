package no.uib.triogen.cmd.transmission;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.PrintWriter;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import no.uib.triogen.TrioGen;
import static no.uib.triogen.io.Utils.lineSeparator;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.readers.SimpleGzReader;
import no.uib.triogen.io.vcf.VcfIterator;
import no.uib.triogen.io.vcf.VcfLine;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.transmission.extraction.Extractor;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;

/**
 * This class writes the summary for a given score on a set of vcf files.
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

            if (!bean.test) {

                run(bean);

            } else {

                runTest(bean);

            }

        } catch (Throwable e) {

            e.printStackTrace();
        }
    }

    /**
     * Runs the command.
     *
     * @param bean the bean of command line parameters
     */
    private static void run(ExtractTransmissionOptionsBean bean) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        Extractor extractor = new Extractor(bean.vcfFile, childToParentMap);

        try {

            extractor.run(
                    bean.nThreads,
                    bean.timeOut
            );

        } catch (Throwable e) {

            e.printStackTrace();

        }
    }

    /**
     * Runs a test.
     *
     * @param bean the bean of command line parameters
     */
    private static void runTest(ExtractTransmissionOptionsBean bean) {

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        SimpleFileWriter writer = new SimpleFileWriter(bean.destinationFile, true);

        VcfIterator vcfIterator = new VcfIterator(bean.vcfFile);

        long start = Instant.now().getEpochSecond();

        for (int i = 0; i < 1000; i++) {

            VcfLine vcfLine = vcfIterator.next();
            vcfLine.parse();

            StringBuilder sb = new StringBuilder();
            sb.append(vcfLine.getVariantDescription());

            for (String childId : childToParentMap.children) {

                String motherId = childToParentMap.getMother(childId);
                String fatherId = childToParentMap.getFather(childId);

                int genotypeChild = vcfLine.getGenotype(childId);
                int genotypeMother = vcfLine.getGenotype(motherId);
                int genotypeFather = vcfLine.getGenotype(fatherId);

                int total = genotypeChild + genotypeMother + genotypeFather;

                sb.append(Integer.toString(total));

            }
            
            writer.writeLine(sb.toString());
            
        }
        
        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        System.out.println("Read and wrote 1000 lines in " + duration + " seconds.");

        vcfIterator.close();
        writer.close();

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
