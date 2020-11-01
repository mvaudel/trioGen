package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import java.util.HashMap;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * This script generates a target file from a Bolt results file
 *
 * @author Marc Vaudel
 */
public class TriogenToLdhub {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File triogenFolder = new File("/mnt/work/marc/moba/run/triogen/lm_output/gwas/z_placenta_weight");
        File ldhubFolder = new File("/mnt/work/marc/moba/run/ld_hub/z_placenta_weight/input");

        String[] variables = new String[]{"cmf.Bc", "cmf.Bm", "cmf.Bf", "cmf_mt.Bc", "cmf_mt.Bm", "cmf_mt.Bf", "cmf_mt.Bmt"};

        HashMap<String, SimpleFileWriter> writers = new HashMap<>(variables.length);

        for (String variable : variables) {

            File destinationFile = new File(ldhubFolder, variable + ".gz");

            SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true);

            writer.writeLine("rsid\tA1\tA2\tN\tP-value\tBeta\tSE");

            writers.put(variable, writer);

        }

        for (int chromosome = 1; chromosome <= 22; chromosome++) {

            System.out.println(Instant.now() + "    Processing chromosome " + chromosome);

            Instant begin = Instant.now();

            File inputFile = new File(triogenFolder, chromosome + ".z_placenta_weight.gz");

            SimpleFileReader reader = SimpleFileReader.getFileReader(inputFile);

            String line = reader.readLine();

            line = reader.readLine();
            String[] lineSplit = line.split("\t");
            HashMap<String, Integer> betaColumn = new HashMap<>(variables.length);
            HashMap<String, Integer> pColumn = new HashMap<>(variables.length);
            HashMap<String, Integer> seColumn = new HashMap<>(variables.length);

            for (String variable : variables) {

                for (int i = 0; i < lineSplit.length; i++) {

                    String colName = lineSplit[i];

                    if (colName.equals(variable)) {

                        betaColumn.put(variable, i);

                    } else if (colName.equals(variable + ".p")) {

                        pColumn.put(variable, i);

                    } else if (colName.equals(variable + ".se")) {

                        seColumn.put(variable, i);

                    }
                }
                
                if (!betaColumn.containsKey(variable)) {
                    
                    System.out.println(line);
                    throw new IllegalArgumentException(variable + " not found in chromosome " + chromosome);
                    
                }
                
                if (!pColumn.containsKey(variable)) {
                    
                    System.out.println(line);
                    throw new IllegalArgumentException(variable + ".p not found in chromosome " + chromosome);
                    
                }
                
                if (!seColumn.containsKey(variable)) {
                    
                    System.out.println(line);
                    throw new IllegalArgumentException(variable + ".se not found in chromosome " + chromosome);
                    
                }
                
            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split("\t");

                String rsid = lineSplit[3];
                String a1 = lineSplit[4];
                String a2 = lineSplit[5];
                String n = lineSplit[7];

                for (String variable : variables) {

                    double beta = Double.parseDouble(lineSplit[betaColumn.get(variable)]);
                    double p = Double.parseDouble(lineSplit[pColumn.get(variable)]);
                    double se = Double.parseDouble(lineSplit[seColumn.get(variable)]);

                    if (!Double.isNaN(beta) && !Double.isInfinite(beta)
                            && !Double.isNaN(p) && !Double.isInfinite(p)
                            && !Double.isNaN(se) && !Double.isInfinite(se)) {

                        String exportLine = String.join("\t",
                                rsid,
                                a1,
                                a2,
                                n,
                                Double.toString(p),
                                Double.toString(beta),
                                Double.toString(se)
                        );

                        SimpleFileWriter writer = writers.get(variable);

                        writer.writeLine(exportLine);

                    }
                }
            }

            Instant end = Instant.now();

            long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println(Instant.now() + "    Chromosome " + chromosome + " finished (" + durationSeconds + " s)");

        }

        for (SimpleFileWriter writer : writers.values()) {

            writer.close();

        }
    }
}
