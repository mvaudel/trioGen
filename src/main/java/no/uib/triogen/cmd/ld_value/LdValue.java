package no.uib.triogen.cmd.ld_value;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import static no.uib.triogen.io.IoUtils.LINE_SEPARATOR;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.processing.ld.LdMatrixComputer;
import no.uib.triogen.utils.Utils;
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
        
        Instant begin = Instant.now();
        
        System.out.println("Initiating LD files.");
        
        boolean wildCard = bean.ldMatrixFilePath.contains(Utils.CONTIG_WILDCARD);
        
        VariantList variantList = bean.variantFile == null ? null : VariantList.getVariantList(bean.variantFile);
        
        HashMap<String, LdMatrixReader> ldMatrixReaderMap = new HashMap<>();
        
        for (String contig : variantList.chromosome) {
            
            if (!ldMatrixReaderMap.containsKey(contig)) {
                
                String ldMatrixFilePath = wildCard ? bean.ldMatrixFilePath.replace(Utils.CONTIG_WILDCARD, contig) : bean.ldMatrixFilePath;
                
                File ldMatrixFile = new File(ldMatrixFilePath);
                
                if (!ldMatrixFile.exists()) {
                    
                    throw new FileNotFoundException("LD matrix file not found: " + ldMatrixFilePath + ".");
                    
                }
                
                LdMatrixReader ldMatrixReader = new LdMatrixReader(ldMatrixFile);                
                
                ldMatrixReaderMap.put(contig, ldMatrixReader);
                
            }
            
        }
        
        Instant end = Instant.now();
        
        long timeInSec = end.getEpochSecond() - begin.getEpochSecond();
        
        System.out.println("Initiating LD files finished (" + timeInSec + " s)");
        
        System.out.println("Getting LD for " + variantList.variantId.length + " variants.");
        
        begin = Instant.now();
        
        Map<String, HashMap<String, Double>> results = IntStream.range(0, variantList.variantId.length)
                .parallel()
                .mapToObj(i -> i)
                .collect(
                        Collectors.toMap(
                                i -> variantList.variantId[i],
                                i -> ldMatrixReaderMap.get(variantList.chromosome[i]).getR2(LINE_SEPARATOR)
                        )
                );
        
        end = Instant.now();
        
        timeInSec = end.getEpochSecond() - begin.getEpochSecond();
        
        System.out.println("Getting LD finished (" + timeInSec + " s)");
        
        System.out.println("Writing results to " + bean.destinationFile + ".");
        
        begin = Instant.now();
        
        JSONObject resultJson = new JSONObject();
        
        for (Entry<String, HashMap<String, Double>> entry : results.entrySet()) {
            
            HashMap<String, Double> ldMap = entry.getValue();
            
            if (ldMap != null) {
                
                String variantId = entry.getKey();
                JSONObject ldMapJson = new JSONObject(ldMap);
                
                resultJson.put(variantId, ldMapJson);
                
            }
        }
        
        FileWriter writer = new FileWriter(bean.destinationFile);
        writer.write(resultJson.toString());
        
        end = Instant.now();
        
        timeInSec = end.getEpochSecond() - begin.getEpochSecond();
        
        System.out.println("Writing results finished (" + timeInSec + " s)");
        
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
