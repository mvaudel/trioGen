package no.uib.triogen.scripts_marc;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts lines corresponding to a set of variants.
 *
 * @author Marc Vaudel
 */
public class ExtractVariants {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File refFile = new File("C:\\Github\\helgeland2020\\resources\\snp-to-locusname-150920.txt");

            HashSet<String> variantIds = new HashSet<>();

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(refFile)) {

                String line = reader.readLine();

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    variantIds.add(lineSplit[0]);

                }
            }

            File proxyFolder = new File("C:\\Github\\helgeland2020\\resources\\warrington");

            HashMap<String, HashMap<String, Double>> proxyR2Map = new HashMap<>();
            HashMap<String, HashMap<String, String>> proxyAlleleMap = new HashMap<>();

            for (File proxyFile : proxyFolder.listFiles()) {

                String fileName = proxyFile.getName();

                if (fileName.startsWith("proxy_")) {

                    String variantId = fileName.substring(6, fileName.length() - 4);

                    try (SimpleFileReader reader = SimpleFileReader.getFileReader(proxyFile)) {

                        String line = reader.readLine();

                        while ((line = reader.readLine()) != null) {

                            reader.readLine();

                            String[] lineSplit = line.split("\t");

                            double r2 = Double.parseDouble(lineSplit[6]);

                            if (r2 >= 0.2) {

                                String proxyId = lineSplit[0];

                                if (proxyId.startsWith("rs")) {

                                    String alleles = lineSplit[7];

                                    HashMap<String, Double> proxyR2SubMap = proxyR2Map.get(proxyId);
                                    HashMap<String, String> proxyAlleleSubMap = proxyAlleleMap.get(proxyId);

                                    if (proxyR2SubMap == null) {

                                        proxyR2SubMap = new HashMap<>(2);
                                        proxyR2Map.put(proxyId, proxyR2SubMap);

                                        proxyAlleleSubMap = new HashMap<>(2);
                                        proxyAlleleMap.put(proxyId, proxyAlleleSubMap);

                                    }

                                    proxyR2SubMap.put(variantId, r2);
                                    proxyAlleleSubMap.put(variantId, alleles);

                                }
                            }
                        }
                    }
                }
            }

            HashMap<String, Double> r2Map = new HashMap<>();
            HashMap<String, String> linesMap = new HashMap<>();

            String headerLine = null;

            File gzFile = new File("C:\\Github\\helgeland2020\\resources\\warrington\\Fetal_BW_European_meta.NG2019.txt.gz");

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(gzFile)) {

                String line = reader.readLine();

                headerLine = line;

                while ((line = reader.readLine()) != null) {

                    reader.readLine();

                    String[] lineSplit = line.split(" ");

                    String proxyId = lineSplit[lineSplit.length - 1];

                    if (proxyId.startsWith("rs") && variantIds.contains(proxyId)) {

                        String newLine = String.join(" ",
                                line,
                                proxyId,
                                "1.0",
                                "-"
                        );

                        linesMap.put(proxyId, newLine);

                        r2Map.put(proxyId, 1.0);

                    }

                    if (proxyAlleleMap.containsKey(proxyId)) {

                        for (Entry<String, Double> r2Entry : proxyR2Map.get(proxyId).entrySet()) {

                            String variantId = r2Entry.getKey();
                            double r2 = r2Entry.getValue();

                            Double foundR2 = r2Map.get(variantId);

                            if (foundR2 == null || foundR2 < r2) {

                                String alleles = proxyAlleleMap.get(proxyId).get(variantId);

                                String newLine = String.join(" ",
                                        line,
                                        variantId,
                                        Double.toString(r2),
                                        alleles
                                );

                                linesMap.put(variantId, newLine);

                                r2Map.put(variantId, r2);

                            }
                        }
                    }
                }
            }

            File destinationFile = new File("C:\\Github\\helgeland2020\\resources\\warrington\\Fetal_BW_European_meta.NG2019_H2020.gz");

            try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

                writer.writeLine(headerLine + " snp_h2020 r2 alleles");

                File missingFile = new File("C:\\Github\\helgeland2020\\resources\\warrington\\missing");

                try (SimpleFileWriter writerMissing = new SimpleFileWriter(missingFile, false)) {

                    for (String variantId : variantIds) {

                        String line = linesMap.get(variantId);

                        if (line != null) {

                            writer.writeLine(line);

                        } else {

                            writerMissing.writeLine(variantId);

                        }
                    }
                }
            }

        } catch (Exception e) {

            e.printStackTrace();

        }
    }
}
