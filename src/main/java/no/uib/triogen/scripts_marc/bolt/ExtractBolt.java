package no.uib.triogen.scripts_marc.bolt;

import java.io.File;
import java.time.Instant;
import java.util.HashMap;
import java.util.stream.IntStream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * Extracts variants from bolt files around a set of targets.
 *
 * @author Marc Vaudel
 */
public class ExtractBolt {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // Get targets
        System.out.println(Instant.now() + "    Loading variants.");

        Instant begin = Instant.now();

        File targetFile = new File("/mnt/work/marc/moba/helgeland2020/scripts_figures/multi-signals/gwas_results/targets");
        VariantList variantList = VariantList.getVariantList(targetFile);

        int maxDistance = 500000;

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Variants lodaed, " + variantList.variantId.length + " variants to map (" + durationSeconds + " s)");

        // Process BOLT file
        String boltFilePath = "/mnt/work2/helgeland/helgeland2020/results/child-additive/z_bmi{ageI}-stats-bgen.gz";
        String resultsFilePath = "/mnt/work/marc/moba/helgeland2020/scripts_figures/multi-signals/gwas_results/{variantId}_zBMI{ageI}.gz";

        IntStream.range(0, 12)
                .parallel()
                .forEach(
                        ageI -> processAgeI(ageI, variantList, maxDistance, boltFilePath, resultsFilePath)
                );

    }

    private static void processAgeI(
            int ageI,
            VariantList variantList,
            int maxDistance,
            String filePathPattern,
            String resultsFilePathPattern
    ) {

        System.out.println(Instant.now() + "    Processing results for z_bmi" + ageI + ".");

        Instant begin = Instant.now();

        String filePath = filePathPattern.replace("{ageI}", Integer.toString(ageI));

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(new File(filePath))) {

            String line = reader.readLine();

            HashMap<String, SimpleFileWriter> writerMap = new HashMap<>(variantList.variantId.length);

            for (String variantId : variantList.variantId) {

                String resultsFilePath = resultsFilePathPattern
                        .replace("{ageI}", Integer.toString(ageI))
                        .replace("{variantId}", variantId);

                SimpleFileWriter writer = new SimpleFileWriter(new File(resultsFilePath), true);

                writer.writeLine(line);

                writerMap.put(variantId, writer);

            }

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                String chr = lineSplit[1];
                int pos = Integer.parseInt(lineSplit[2]);

                for (int variantI = 0; variantI < variantList.variantId.length; variantI++) {

                    if (chr.equals(variantList.chromosome[variantI]) && Math.abs(pos - variantList.position[variantI]) <= maxDistance) {

                        SimpleFileWriter writer = writerMap.get(variantList.variantId[variantI]);

                        writer.writeLine(line);

                    }
                }
            }

            for (SimpleFileWriter writer : writerMap.values()) {

                writer.close();

            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    " + ageI + " finished (" + durationSeconds + " s)");

    }
}
