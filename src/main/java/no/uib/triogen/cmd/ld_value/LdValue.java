package no.uib.triogen.cmd.ld_value;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.TrioGen;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.utils.Utils;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * Computes LD between variants and saves the results in a matrix.
 *
 * @author Marc Vaudel
 */
public class LdValue {

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
            LdValueOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            LdValueOptionsBean bean = new LdValueOptionsBean(commandLine);

            run(
                    bean
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
     * @throws IOException Exception thrown if an error occurs while reading the
     * ld file.
     */
    private static void run(
            LdValueOptionsBean bean
    ) throws IOException {

        VariantList variantList = VariantList.getVariantList(bean.variantFile);

        Instant begin = Instant.now();

        System.out.println("Initiating LD files.");

        boolean wildCard = bean.ldMatrixFilePath.contains(Utils.CONTIG_WILDCARD);

        ConcurrentHashMap<String, LdMatrixReader> ldMatrixReaderMap = new ConcurrentHashMap<>();

        HashSet<String> contigs = Arrays.stream(variantList.chromosome)
                .collect(
                        Collectors.toCollection(HashSet::new)
                );

        if (wildCard) {

            contigs.parallelStream()
                    .forEach(
                            contig -> {

                                String ldMatrixFilePath = bean.ldMatrixFilePath.replace(Utils.CONTIG_WILDCARD, contig);

                                File ldMatrixFile = new File(ldMatrixFilePath);

                                if (!ldMatrixFile.exists()) {

                                    throw new IllegalArgumentException("LD matrix file not found: " + ldMatrixFilePath + ".");

                                }

                                try {

                                    LdMatrixReader ldMatrixReader = new LdMatrixReader(ldMatrixFile);

                                    ldMatrixReaderMap.put(contig, ldMatrixReader);

                                } catch (Exception e) {

                                    throw new RuntimeException(e);

                                }
                            }
                    );
        } else {

            LdMatrixReader ldMatrixReader = new LdMatrixReader(new File(bean.ldMatrixFilePath));

            contigs.forEach(
                    contig -> ldMatrixReaderMap.put(contig, ldMatrixReader)
            );
        }

        Instant end = Instant.now();

        long timeInSec = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println("Initiating LD files finished (" + ldMatrixReaderMap.size() + " files in " + timeInSec + " s)");

        System.out.println("Getting LD for " + variantList.variantId.length + " variants.");

        begin = Instant.now();

        ConcurrentHashMap<String, ArrayList<R2>> results = new ConcurrentHashMap<>(variantList.variantId.length);

        IntStream.range(0, variantList.variantId.length)
                .parallel()
                .forEach(
                        i -> {
                            String variantId = variantList.variantId[i];
                            String contig = variantList.chromosome[i];
                            LdMatrixReader ldMatrixReader = ldMatrixReaderMap.get(contig);
                            ArrayList<R2> result = ldMatrixReader.getR2(variantId);

                            if (result != null) {

                                for (R2 r2 : result) {

                                    String variantB = ldMatrixReader.getId(r2.variantB);
                                    r2.setVariantBId(variantB);

                                    String rsidB = ldMatrixReader.getRsId(variantB);
                                    r2.setVariantBRsid(rsidB);

                                }

                                results.put(variantId, result);

                            }
                        }
                );

        end = Instant.now();

        timeInSec = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println("Getting LD finished (" + timeInSec + " s)");

        System.out.println("Writing results to " + bean.destinationFile + ".");

        begin = Instant.now();

        JSONObject resultJson = new JSONObject();

        for (String variantId : variantList.variantId) {

            ArrayList<R2> r2Array = results.get(variantId);

            if (r2Array != null) {

                JSONArray array = new JSONArray();

                for (R2 r2 : r2Array) {

                    JSONObject r2Json = new JSONObject();

                    r2Json.put("variantA", variantId);
                    r2Json.put("alleleA", r2.alleleA);
                    r2Json.put("variantB", r2.getVariantBId());
                    r2Json.put("rsidB", r2.getVariantBRsid());
                    r2Json.put("alleleB", r2.alleleB);
                    r2Json.put("r2", r2.r2Value);

                    array.put(r2Json);

                }

                resultJson.put(variantId, array);

            } else {

                resultJson.put(variantId, "Not Found");

            }
        }

        try (SimpleFileWriter writer = new SimpleFileWriter(bean.destinationFile, true)) {

            writer.writeLine(resultJson.toString());

        }

        end = Instant.now();

        timeInSec = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println("Writing results finished (" + timeInSec + " s)");

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
            lPrintWriter.print("    Linkage Disequilibrium Value  " + LINE_SEPARATOR);
            lPrintWriter.print("==================================" + LINE_SEPARATOR);
            lPrintWriter.print(LINE_SEPARATOR
                    + "The Linkage Disequilibrium Value command line returns the variants in ld with a given set of variants." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + LINE_SEPARATOR
                    + LINE_SEPARATOR
                    + "----------------------"
                    + LINE_SEPARATOR
                    + "OPTIONS"
                    + LINE_SEPARATOR
                    + "----------------------" + LINE_SEPARATOR
                    + LINE_SEPARATOR);
            lPrintWriter.print(LdValueOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
