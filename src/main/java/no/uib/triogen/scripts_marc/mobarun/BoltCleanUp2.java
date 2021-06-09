package no.uib.triogen.scripts_marc.mobarun;

import java.io.File;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.Instant;
import java.util.ArrayList;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;

/**
 * This script Cleans the organization of BOLT-LMM files on MoBa.
 *
 * @author Marc Vaudel
 */
public class BoltCleanUp2 {

    public final static double MAF_THRESHOLD = 0.001;

    private static int trimmingProgress = 0;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            // Iterate all phenos check which analyses are finished, which are incomplete or failed
            File boltResultsFolder = new File("/mnt/work/marc/moba/run/bolt/bolt_output");

            String[] genos = new String[]{"child", "mother", "father"};

            String[] boltFileTemplates = new String[]{"{geno}Geno_{pheno}-stats-bgen.gz", "{geno}Geno_{pheno}-stats-bgen-chrX.gz"};
            String[] boltTabFileTemplates = new String[]{"{geno}Geno_{pheno}-stats.tab", "{geno}Geno_{pheno}-stats-chrX.tab"};
            String[] logFileTemplates = new String[]{"{geno}Geno_{pheno}-runlog.log", "{geno}Geno_{pheno}-runlog-chrX.log"};
            String[] rareFileTemplates = new String[]{"{geno}Geno_{pheno}_under_maf_" + MAF_THRESHOLD + ".gz", "{geno}Geno_{pheno}_under_maf_" + MAF_THRESHOLD + "-chrX.gz"};
            String[] commonFileTemplates = new String[]{"{geno}Geno_{pheno}_over_maf_" + MAF_THRESHOLD + ".gz", "{geno}Geno_{pheno}_over_maf_" + MAF_THRESHOLD + "-chrX.gz"};
            String[] chrLabels = new String[]{"1:22", "X"};

            for (File phenoFolder : boltResultsFolder.listFiles()) {

                if (phenoFolder.isDirectory()) {

                    String pheno = phenoFolder.getName();

                    System.out.println(Instant.now() + "    Cleaning up old files for " + pheno + ".");

                    File prunedFolder = new File(phenoFolder, "pruned");

                    if (prunedFolder.exists()) {

                        deleteDirContent(prunedFolder);

                    }

                    System.out.println(Instant.now() + "    Inspecting files for " + pheno + ".");

                    Instant begin = Instant.now();

                    StringBuilder report = new StringBuilder();

                    for (String geno : genos) {

                        for (int chrI = 0; chrI < 2; chrI++) {

                            String rareFileName = rareFileTemplates[chrI]
                                    .replace("{geno}", geno)
                                    .replace("{pheno}", pheno);
                            String commonFileName = commonFileTemplates[chrI]
                                    .replace("{geno}", geno)
                                    .replace("{pheno}", pheno);

                            File rareFile = new File(phenoFolder, rareFileName);
                            File rareFileCopy = new File(phenoFolder, "temp1");

                            File commonFile = new File(phenoFolder, commonFileName);
                            File commonFileCopy = new File(phenoFolder, "temp2");

                            if (rareFile.exists()) {

                                Files.move(rareFile.toPath(), rareFileCopy.toPath(), StandardCopyOption.ATOMIC_MOVE);

                            }

                            if (commonFile.exists()) {

                                Files.move(commonFile.toPath(), commonFileCopy.toPath(), StandardCopyOption.ATOMIC_MOVE);

                            }

                            if (rareFileCopy.exists()) {

                                Files.move(rareFileCopy.toPath(), commonFile.toPath(), StandardCopyOption.ATOMIC_MOVE);

                            }

                            if (commonFileCopy.exists()) {

                                Files.move(commonFileCopy.toPath(), rareFile.toPath(), StandardCopyOption.ATOMIC_MOVE);

                            }
                        }
                    }

                    Instant end = Instant.now();
                    long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                    System.out.println(Instant.now() + "    " + pheno + ": " + report + " (" + durationSeconds + " s)");

                }
            }
        } catch (Throwable e) {

            e.printStackTrace();

        }
    }

    private static void deleteDir(File dir) {

        if (dir.isDirectory()) {

            deleteDirContent(dir);

        }

        dir.delete();

    }

    private static void deleteDirContent(File dir) {

        if (dir.isDirectory()) {

            for (File file : dir.listFiles()) {

                deleteDir(file);

            }
        }
    }
}
