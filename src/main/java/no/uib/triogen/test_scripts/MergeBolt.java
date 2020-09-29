/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uib.triogen.test_scripts;

import java.io.File;
import java.time.Instant;
import java.util.HashSet;
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

        HashSet<String> targets = new HashSet<>();

        String[] targetPaths = new String[]{
            "/mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets_h2020",
            "/mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets"
        };

        for (String targetPath : targetPaths) {

            File targetsFile = new File(targetPath);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(targetsFile)) {

                String line = reader.readLine();

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    targets.add(lineSplit[0]);

                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Importing targets finished (" + targets.size() + " imported in " + durationSeconds + " s)");

        // Iterate bolt results and save target variants
        File destinationFile = new File("/mnt/work/marc/moba/H2020/docs/results/parent_gwas.gz");

        boolean header = false;

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            String[] phenos = new String[]{
                "z_mother_median_height", "z_mother_weight_beginning", "z_mother_bmi_beginning", "z_father_height", "z_father_weight", "z_father_bmi"
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

                                    if (targets.contains(rsId)) {

                                        line = String.join("\t", pheno, geno, line);

                                        writer.writeLine(line);

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
    }

}
