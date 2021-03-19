package no.uib.triogen.scripts_marc.mobarun;

import java.io.File;
import java.io.IOException;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.Instant;
import java.util.Arrays;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;

/**
 * This script Cleans the organization of BOLT-LMM files on MoBa.
 *
 * @author Marc Vaudel
 */
public class BoltCleanUp {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // Check that every file is in its own folder
        File boltResultsFolder = new File("/mnt/work/marc/moba/run/bolt/bolt_output");
        String[] genos = new String[]{"child", "mother", "father"};
        String[] fileTemplates = new String[]{
            "{geno}_{pheno}-runlog.log",
            "{geno}_{pheno}-runlog-chrX.log",
            "{geno}_{pheno}-stats-bgen.gz",
            "{geno}_{pheno}-stats-bgen-chrX.gz",
            "{geno}_{pheno}-stats.tab",
            "{geno}_{pheno}-stats-chrX.tab"
        };

        for (File file : boltResultsFolder.listFiles()) {

            if (file.exists()) {

                String fileName = file.getName();

                if (fileName.endsWith("stats-bgen.gz")) {

                    String pheno = fileName
                            .substring(fileName.indexOf("_") + 1, fileName.length() - 14);

                    for (String geno : genos) {

                        System.out.println(Instant.now() + "    Processing " + pheno + ".");

                        Instant begin = Instant.now();

                        File phenoFolder = new File(boltResultsFolder, pheno);

                        if (!phenoFolder.exists()) {

                            phenoFolder.mkdir();

                        }

                        Arrays.stream(fileTemplates)
                                .parallel()
                                .forEach(
                                        fileTemplate -> processBoltFile(fileName, geno, pheno, phenoFolder)
                                );

                        Instant end = Instant.now();

                        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                        System.out.println(Instant.now() + "    Processing " + pheno + " finished (" + durationSeconds + " s)");

                    }
                }
            }
        }
    }

    private static void processBoltFile(
            String fileTemplate,
            String geno,
            String pheno,
            File phenoFolder
    ) {

        try {

            String boltFilePath = fileTemplate
                    .replace("{geno}", geno)
                    .replace("{pheno}", pheno);

            File boltFile = new File(boltFilePath);
            File destinationFile = new File(phenoFolder, boltFile.getName());

            if (boltFile.exists()) {

                Files.move(
                        boltFile.toPath(), 
                        destinationFile.toPath(), 
                        StandardCopyOption.ATOMIC_MOVE, 
                        StandardCopyOption.REPLACE_EXISTING
                );

            }

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }
}
