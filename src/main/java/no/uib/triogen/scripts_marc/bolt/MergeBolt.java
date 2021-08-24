/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uib.triogen.scripts_marc.bolt;

import java.io.File;
import java.time.Instant;
import java.util.HashMap;
import java.util.Map.Entry;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Script to merge BOLT results from different phenotypes.
 *
 * @author mvaudel
 */
public class MergeBolt {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // get targets
        System.out.println(Instant.now() + "    Importing targets");

        Instant begin = Instant.now();

        HashMap<String, Boolean> targets = new HashMap<>();

        String[] targetPaths = new String[]{
            "/mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets_egg_bw",
            "/mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets_yengo_bmi"
        };

        for (String targetPath : targetPaths) {

            File targetsFile = new File(targetPath);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(targetsFile)) {

                String line = reader.readLine();

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    targets.put(lineSplit[0], false);

                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Importing targets finished (" + targets.size() + " imported in " + durationSeconds + " s)");

        // Iterate bolt results and save target variants
        File destinationFile = new File("/mnt/work/marc/moba/H2020/docs/results/yengo_bw.gz");

        boolean header = false;

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            String[] phenos = new String[]{
                "z_bmi0", 
                "z_bmi1", 
                "z_bmi2", 
                "z_bmi3", 
                "z_bmi4", 
                "z_bmi5", 
                "z_bmi6", 
                "z_bmi7", 
                "z_bmi8", 
                "z_bmi9", 
                "z_bmi10", 
                "z_bmi11"
            };

            String[] genomes = new String[]{
                "child", "mother", "father"
            };

            String[] suffixes = new String[]{
                "", "-chrX"
            };

            for (int i = 0; i < phenos.length; i++) {

                String pheno = phenos[i];

                System.out.println(Instant.now() + "    Processing " + pheno + " (" + (i + 1) + " of " + phenos.length + ")");

                begin = Instant.now();

                for (String suffix : suffixes) {

                    for (String geno : genomes) {

                        String boltFilePath = "/mnt/work/marc/moba/run/bolt/bolt_output/" + geno + "Geno_" + pheno + "-stats-bgen" + suffix + ".gz";

                        File boltFile = new File(boltFilePath);

                        if (boltFile.exists()) {

                            try (SimpleFileReader reader = SimpleFileReader.getFileReader(boltFile)) {

                                String line = reader.readLine();

                                if (!header) {

                                    String headerLine = String.join("\t", "PHENO", "GENO", line);

                                    writer.writeLine(headerLine);

                                    header = true;

                                }

                                while ((line = reader.readLine()) != null) {

                                    int separatorIndex = line.indexOf('\t');

                                    String rsId = line.substring(0, separatorIndex);

                                    if (targets.containsKey(rsId)) {

                                        line = String.join("\t", pheno, geno, line);

                                        writer.writeLine(line);

                                        targets.put(rsId, true);

                                    }
                                }
                            }
                        }
                    }
                }

                end = Instant.now();

                durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                System.out.println(Instant.now() + "    " + pheno + " finished (" + durationSeconds + " s)");

            }
        }

        String[] missing = targets.entrySet()
                .stream()
                .filter(
                        entry -> entry.getValue() == false
                )
                .map(
                        Entry::getKey
                )
                .toArray(String[]::new);

        if (missing.length > 0) {
            
            System.out.println(missing.length + " variants not found in MoBa.");

            File missingFile = new File("/mnt/work/marc/moba/H2020/docs/results/yengo_bw_missing.gz");

            try (SimpleFileWriter missingWriter = new SimpleFileWriter(missingFile, true)) {

                for (String key : missing) {

                    missingWriter.writeLine(key);

                }
            }
        } else {
            
            System.out.println("All variants found in MoBa.");
            
        }
    }
}
