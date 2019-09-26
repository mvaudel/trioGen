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
import no.uib.triogen.io.flat.readers.SimpleGzReader;
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

        VCFFileReader vcfFileReader = new VCFFileReader(bean.vcfFile);
        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(bean.trioFile);

        Extractor extractor = new Extractor(vcfFileReader, childToParentMap);

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

        String childId = childToParentMap.children.first();
        String motherId = childToParentMap.getMother(childId);
        String fatherId = childToParentMap.getFather(childId);
        
        SimpleGzReader reader = new SimpleGzReader(bean.vcfFile);

        VCFFileReader vcfFileReader = new VCFFileReader(bean.vcfFile);
        CloseableIterator<VariantContext> iterator = vcfFileReader.iterator();
        
        VariantContext variantContext;

        while ((variantContext = iterator.next()) != null) {

            ArrayList<String> alleles = variantContext.getAlleles().stream()
                    .map(
                            allele -> allele.getBaseString()
                    )
                    .collect(
                            Collectors.toCollection(ArrayList::new)
                    );

            if (!alleles.get(0).equals(alleles.get(1))) {

                Genotype childGenotype = variantContext.getGenotype(childId);
                List<Allele> childAlleles = childGenotype.getAlleles();
                String childAllelesString = childAlleles.stream()
                        .map(
                                allele -> allele.getBaseString()
                        )
                        .collect(
                                Collectors.joining(",")
                        );

                Genotype motherGenotype = variantContext.getGenotype(motherId);
                List<Allele> motherAlleles = motherGenotype.getAlleles();
                String motherAllelesString = motherAlleles.stream()
                        .map(
                                allele -> allele.getBaseString()
                        )
                        .collect(
                                Collectors.joining(",")
                        );

                Genotype fatherGenotype = variantContext.getGenotype(fatherId);
                List<Allele> fatherAlleles = fatherGenotype.getAlleles();
                String fatherAllelesString = fatherAlleles.stream()
                        .map(
                                allele -> allele.getBaseString()
                        )
                        .collect(
                                Collectors.joining(",")
                        );

                if (!motherAllelesString.equals(childAllelesString)) {

                    System.out.println("Id: " + variantContext.getID());
                    System.out.println("Ref: " + alleles.get(0) + " Alt: " + alleles.get(1));
                    System.out.println("Mother genotype: " + motherAllelesString);
                    System.out.println("Child genotype: " + childAllelesString);
                    System.out.println("Father genotype: " + fatherAllelesString);

                    break;
                }
            }
        }

        long start = Instant.now().getEpochSecond();

        int nVariants = 0;

        String line = reader.readLine();
        while ((line = reader.readLine()) != null && ++nVariants < 1000) {
        }

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        System.out.println("Read " + nVariants + " lines in " + duration + " seconds.");

        start = Instant.now().getEpochSecond();

        nVariants = 0;

        while ((variantContext = iterator.next()) != null && ++nVariants < 1000) {
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        System.out.println("Read " + nVariants + " variants in " + duration + " seconds.");

        start = Instant.now().getEpochSecond();

        nVariants = 0;

        while ((variantContext = iterator.next()) != null && ++nVariants < 1000) {

            for (String tempChildId : childToParentMap.children) {

                String tempMotherId = childToParentMap.getMother(tempChildId);
                String tempFatherId = childToParentMap.getFather(tempChildId);

                variantContext.getGenotype(tempChildId);
                variantContext.getGenotype(tempMotherId);
                variantContext.getGenotype(tempFatherId);

            }
        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        System.out.println("Processed " + nVariants + " variants in " + duration + " seconds.");

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
