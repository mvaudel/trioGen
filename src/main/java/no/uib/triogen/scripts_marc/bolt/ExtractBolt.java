package no.uib.triogen.scripts_marc.bolt;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
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
     * The maximal distance allowed around a hit (inclusive).
     */
    private static final int MAX_DISTANCE = 0;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // Get targets
        System.out.println(Instant.now() + "    Loading variants.");

        Instant begin = Instant.now();

        File targetFile = new File("/mnt/work/marc/moba/puberty_2020/resources/targets/targets_with_proxies_18.11.21");
        VariantList variantList = VariantList.getVariantList(targetFile);

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Variants lodaed, " + variantList.variantId.length + " variants to map (" + durationSeconds + " s)");

        // Process BOLT file
        String boltFilePath = "/mnt/work/marc/moba/run/bolt/bolt_output/childGeno_{pheno}{ageI}{mafSuffix}{chrSuffix}.gz";
        String resultsFilePath = "/mnt/work/marc/moba/puberty_2020/resources/bolt/{pheno}{ageI}.gz";

        String[] phenos = new String[]{"z_bmi", "z_length", "z_weight"};
        String[] mafSuffixes = new String[]{"_over_maf_0.001", "_under_maf_0.001"};
        String[] chrSuffixes = new String[]{"", "-chrX"};

        Arrays.stream(phenos)
                .parallel()
                .forEach(
                        pheno -> IntStream.range(0, 12)
                                .parallel()
                                .forEach(
                                        ageI -> matchById(
                                                pheno,
                                                ageI,
                                                mafSuffixes,
                                                chrSuffixes,
                                                variantList,
                                                boltFilePath,
                                                resultsFilePath
                                        )
                                )
                );

    }

    private static void matchByDistance(
            String pheno,
            int ageI,
            String mafSuffix,
            String chrSuffix,
            VariantList variantList,
            String filePathPattern,
            String resultsFilePathPattern
    ) {

        System.out.println(Instant.now() + "    Processing results for " + pheno + ageI + ".");

        Instant begin = Instant.now();

        String filePath = filePathPattern.replace("{pheno}", pheno);
        filePath = filePath.replace("{ageI}", Integer.toString(ageI));
        filePath = filePath.replace("{mafSuffix}", mafSuffix);
        filePath = filePath.replace("{chrSuffix}", chrSuffix);

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

                    if (chr.equals(variantList.chromosome[variantI]) && Math.abs(pos - variantList.position[variantI]) <= MAX_DISTANCE) {

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

    private static void matchById(
            String pheno,
            int ageI,
            String[] mafSuffixes,
            String[] chrSuffixes,
            VariantList variantList,
            String filePathPattern,
            String resultsFilePathPattern
    ) {

        System.out.println(Instant.now() + "    Processing results for " + pheno + ageI + ".");

        Instant begin = Instant.now();

        String destinationFilePath = resultsFilePathPattern.replace("{pheno}", pheno);
        destinationFilePath = destinationFilePath.replace("{ageI}", Integer.toString(ageI));

        try (SimpleFileWriter writer = new SimpleFileWriter(new File(destinationFilePath), true)) {

            boolean header = false;

            for (String mafSuffix : mafSuffixes) {

                for (String chrSuffix : chrSuffixes) {

                    String boltFilePath = filePathPattern.replace("{pheno}", pheno);
                    boltFilePath = boltFilePath.replace("{ageI}", Integer.toString(ageI));
                    boltFilePath = boltFilePath.replace("{mafSuffix}", mafSuffix);
                    boltFilePath = boltFilePath.replace("{chrSuffix}", chrSuffix);

                    File boltFile = new File(boltFilePath);

                    if (boltFile.exists()) {

                        try (SimpleFileReader reader = SimpleFileReader.getFileReader(boltFile)) {

                            String line = reader.readLine();

                            if (!header) {

                                writer.writeLine(line);

                                header = false;

                            }

                            while ((line = reader.readLine()) != null) {

                                String[] lineSplit = line.split("\t");

                                String snp = lineSplit[0];

                                for (int variantI = 0; variantI < variantList.variantId.length; variantI++) {

                                    if (snp.equals(variantList.variantId[variantI])) {

                                        writer.writeLine(line);

                                    }
                                }
                            }
                        }
                    } else {

                        System.out.println(Instant.now() + "    Bolt file '" + boltFilePath + "' not found");

                    }
                }
            }

            Instant end = Instant.now();

            long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println(Instant.now() + "    " + ageI + " finished (" + durationSeconds + " s)");

        }
    }
}
