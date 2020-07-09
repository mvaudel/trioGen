package no.uib.triogen.test_scripts;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.ld.LdMatrixReader;

/**
 * Extracts the top snps that are not in ld from bolt-lmm results.
 *
 * @author Marc Vaudel
 */
public class PruneBolt {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            double pValueThreshold = 1e-6;
            double ldThreshold = 0.05;

            String[] roles = new String[]{"child", "mother", "father"};
            String[] phenotypes = new String[]{"z_umbilical_cord_length"};

            File ldFolder = new File("/mnt/work/marc/moba/run/triogen/ld");

            File resultsFolder = new File("/mnt/work/marc/moba/run/bolt/bolt_output");

            for (String phenoName : phenotypes) {

                for (String role : roles) {

                    System.out.println("Pruning " + phenoName + " " + role + ".");

                    Instant begin = Instant.now();

                    String previousContig = "";
                    TreeMap<Double, ArrayList<String[]>> pValuesMap = new TreeMap<>();

                    File associationResultsFile = new File(resultsFolder, "MobaRun_" + role + "Geno_" + phenoName + "_maf0.005.gz");
                    File prunedFile = new File(resultsFolder, "MobaRun_" + role + "Geno_" + phenoName + "_maf0.005_pruned.gz");

                    try ( SimpleFileReader reader = SimpleFileReader.getFileReader(associationResultsFile)) {

                        String line = reader.readLine();

                        try ( SimpleFileWriter writer = new SimpleFileWriter(prunedFile, true)) {

                            writer.writeLine(line);

                            while ((line = reader.readLine()) != null) {

                                String[] lineSplit = line.split(" ");

                                String contig = lineSplit[1];

                                if (!contig.equals(previousContig)) {

                                    ArrayList<String[]> prunedLines = prune(pValuesMap, previousContig, ldFolder, ldThreshold);

                                    for (String[] lineSplitPruned : prunedLines) {

                                        String linePruned = String.join(" ", lineSplitPruned);

                                        writer.writeLine(linePruned);

                                    }

                                    pValuesMap.clear();
                                    previousContig = contig;

                                }

                                double p = Double.parseDouble(lineSplit[6]);

                                if (p > 0.0 && p <= pValueThreshold) {

                                    ArrayList<String[]> lines = pValuesMap.get(p);

                                    if (lines == null) {

                                        lines = new ArrayList<>(1);
                                        pValuesMap.put(p, lines);

                                    }

                                    lines.add(lineSplit);

                                }
                            }
                        }
                    }
                }
            }

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

    private static ArrayList<String[]> prune(
            TreeMap<Double, ArrayList<String[]>> pValuesMap,
            String chromosome,
            File ldFolder,
            double ldThreshold
    ) throws IOException {

        ArrayList<String[]> result = new ArrayList<>();

        File ldFile = new File(ldFolder, "dos_chr_" + chromosome + ".tld");

        LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);

        HashSet<String> inspectedSnp = new HashSet<>();

        for (Entry<Double, ArrayList<String[]>> entry : pValuesMap.entrySet()) {

            for (String[] lineSplit : entry.getValue()) {

                String variantId = lineSplit[0];

                if (!inspectedSnp.contains(variantId)) {

                    result.add(lineSplit);
                    inspectedSnp.add(variantId);

                }

                HashMap<String, Double> variantLdMap = ldMatrixReader.getR2(variantId);

                for (Entry<String, Double> entry2 : variantLdMap.entrySet()) {

                    if (entry2.getValue() >= ldThreshold) {

                        inspectedSnp.add(entry2.getKey());

                    }

                }
            }
        }

        return result;

    }

}
