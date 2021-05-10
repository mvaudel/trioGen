package no.uib.triogen.scripts_marc.mobarun;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;

/**
 * This script Cleans the organization of BOLT-LMM files on MoBa.
 *
 * @author Marc Vaudel
 */
public class BoltCleanUp {

    public final static double MAF_THRESHOLD = 0.001;

    private static int trimmingProgress = 0;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // Iterate all phenos check which analyses are finished, which are incomplete or failed
        File boltResultsFolder = new File("/mnt/work/marc/moba/run/bolt/bolt_output");
        File incompleteFileList = new File(boltResultsFolder, "incomplete");
        File failedFileList = new File(boltResultsFolder, "failed");

        String[] genos = new String[]{"child", "mother", "father"};

        String[] boltFileTemplates = new String[]{"{geno}Geno_{pheno}-stats-bgen.gz", "{geno}Geno_{pheno}-stats-bgen-chrX.gz"};
        String[] boltTabFileTemplates = new String[]{"{geno}Geno_{pheno}-stats.tab", "{geno}Geno_{pheno}-stats-chrX.tab"};
        String[] logFileTemplates = new String[]{"{geno}Geno_{pheno}-runlog.log", "{geno}Geno_{pheno}-runlog-chrX.log"};
        String[] rareFileTemplates = new String[]{"{geno}Geno_{pheno}_over_maf_" + MAF_THRESHOLD + ".gz", "{geno}Geno_{pheno}_over_maf_" + MAF_THRESHOLD + "-chrX.gz"};
        String[] commonFileTemplates = new String[]{"{geno}Geno_{pheno}_under_maf_" + MAF_THRESHOLD + ".gz", "{geno}Geno_{pheno}_under_maf_" + MAF_THRESHOLD + "-chrX.gz"};
        String[] chrLabels = new String[]{"1:22", "X"};

        ArrayList<File[]> rawBoltFilesToTrim = new ArrayList<>();

        try (SimpleFileWriter incompleteWriter = new SimpleFileWriter(incompleteFileList, false)) {

            incompleteWriter.writeLine("# Incomplete analyes" + Instant.now());
            incompleteWriter.writeLine("Phenotype", "Individual", "Chromosome");

            try (SimpleFileWriter failedWriter = new SimpleFileWriter(failedFileList, false)) {

                failedWriter.writeLine("# Failed analyses" + Instant.now());
                failedWriter.writeLine("Phenotype", "Individual", "Chromosome");

                for (File phenoFolder : boltResultsFolder.listFiles()) {

                    if (phenoFolder.isDirectory()) {

                        String pheno = phenoFolder.getName();

                        System.out.println(Instant.now() + "    Cleaning up old files for " + pheno + ".");

                        File trimmedFolder = new File(phenoFolder, "maf_0.001");

                        if (trimmedFolder.exists()) {

                            deleteDir(trimmedFolder);

                        }

                        System.out.println(Instant.now() + "    Inspecting files for " + pheno + ".");

                        Instant begin = Instant.now();

                        StringBuilder report = new StringBuilder();

                        for (String geno : genos) {

                            for (int chrI = 0; chrI < 2; chrI++) {

                                String chrLabel = chrLabels[chrI];
                                String logFileName = logFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);
                                String boltFileName = boltFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);
                                String boltTabFileName = boltTabFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);
                                String rareFileName = rareFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);
                                String commonFileName = commonFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);

                                File logFile = new File(phenoFolder, logFileName);
                                File boltFile = new File(phenoFolder, boltFileName);
                                File boltTabFile = new File(phenoFolder, boltTabFileName);

                                if (logFile.exists()) {

                                    if (boltFile.exists() && boltTabFile.exists()) {

                                        File prunedFolder = new File(phenoFolder, "pruned");

                                        if (prunedFolder.exists()) {

                                            deleteDir(prunedFolder);

                                        }

                                        File rareFile = new File(phenoFolder, rareFileName);
                                        File commonFile = new File(phenoFolder, commonFileName);

                                        rawBoltFilesToTrim.add(
                                                new File[]{boltFile, rareFile, commonFile}
                                        );

                                        if (report.length() > 0) {

                                            report.append(", ");

                                        }
                                        report.append(geno).append(" chr").append(chrLabel).append(" to split");

                                    } else {

                                        failedWriter.writeLine(pheno, geno, chrLabel);

                                        if (report.length() > 0) {

                                            report.append(", ");

                                        }
                                        report.append(geno).append(" chr").append(chrLabel).append(" failed");

                                        boltFile.delete();
                                        boltTabFile.delete();

                                    }

                                } else {

                                    incompleteWriter.writeLine(pheno, geno, chrLabel);

                                    if (report.length() > 0) {

                                        report.append(", ");

                                    }
                                    report.append(geno).append(" chr").append(chrLabel).append(" missing");

                                }
                            }
                        }

                        Instant end = Instant.now();
                        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                        System.out.println(Instant.now() + "    " + pheno + ": " + report + " (" + durationSeconds + " s)");

                    }
                }
            }

            // Split bolt results by maf
            rawBoltFilesToTrim.stream()
                    .parallel()
                    .forEach(
                            filePair -> processBoltFile(
                                    filePair[0],
                                    filePair[1],
                                    filePair[2],
                                    rawBoltFilesToTrim.size(),
                                    incompleteWriter
                            )
                    );
        }

    }

    private static void processBoltFile(
            File boltFile,
            File rareFile,
            File commonFile,
            int nFiles,
            SimpleFileWriter incompleteWriter
    ) {

        System.out.println(Instant.now() + "    Splitting " + boltFile + " (" + trimmingProgress++ + " of " + nFiles + ").");

        Instant begin = Instant.now();

        int nVariants = 0;

        boolean success = true;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(boltFile, false)) {

            try (SimpleFileWriter rareWriter = new SimpleFileWriter(rareFile, true)) {

                try (SimpleFileWriter commonWriter = new SimpleFileWriter(commonFile, true)) {

                    String line = reader.readLine();
                    rareWriter.writeLine(line);
                    commonWriter.writeLine(line);

                    while ((line = reader.readLine()) != null) {

                        String[] lineSplit = line.split("\t");

                        if (!lineSplit[11].equals("-nan") && !lineSplit[11].equals("nan")) {

                            double maf = Double.parseDouble(lineSplit[6]);

                            boolean common = false;

                            if (maf >= MAF_THRESHOLD && maf <= 1.0 - MAF_THRESHOLD) {

                                try {

                                    double beta = Double.parseDouble(lineSplit[10]);

                                    if (!Double.isNaN(beta) && !Double.isInfinite(beta)) {

                                        double se = Double.parseDouble(lineSplit[11]);

                                        if (!Double.isNaN(se) && !Double.isInfinite(se)) {

                                            double p = Double.parseDouble(lineSplit[15]);

                                            if (!Double.isNaN(p) && !Double.isInfinite(p)) {

                                                common = true;

                                            }
                                        }
                                    }

                                } catch (Exception e) {

                                    // Parsing error - ignore
                                }

                                if (common) {

                                    commonWriter.writeLine(line);

                                    nVariants++;

                                } else {

                                    rareWriter.writeLine(line);

                                }
                            }
                        }
                    }
                }
            }

        } catch (Exception e) {

            success = false;

        }

        if (success) {

            boltFile.delete();

            Instant end = Instant.now();
            long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println(Instant.now() + "    Splitting " + boltFile + " done, " + nVariants + " common variants (" + durationSeconds + " s)");

        } else {

            boltFile.delete();
            rareFile.delete();
            commonFile.delete();

            incompleteWriter.writeLine("Splitting crash: ", boltFile.getAbsolutePath());

            Instant end = Instant.now();
            long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println(Instant.now() + "    Splitting " + boltFile + " failed (" + durationSeconds + " s)");

        }
    }

    private static void deleteDir(File dir) {

        if (dir.isDirectory()) {

            for (File file : dir.listFiles()) {

                deleteDir(file);

            }
        }

        dir.delete();

    }
}
