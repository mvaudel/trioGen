package no.uib.triogen.scripts_marc.mobarun;

import java.io.File;
import java.io.IOException;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;
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
        File incompleteFileList = new File("/mnt/work/marc/moba/run/bolt/bolt_output/incomplete");
        File failedFileList = new File("/mnt/work/marc/moba/run/bolt/bolt_output/failed");

        String[] genos = new String[]{"child", "mother", "father"};
        String[] boltFileTemplates = new String[]{"{geno}Geno_{pheno}-stats-bgen.gz", "{geno}Geno_{pheno}-stats-bgen-chrX.gz"};
        String[] boltTabFileTemplates = new String[]{"{geno}Geno_{pheno}-stats.tab", "{geno}Geno_{pheno}-stats-chrX.tab"};
        String[] logFileTemplates = new String[]{"{geno}Geno_{pheno}-runlog.log", "{geno}Geno_{pheno}-runlog-chrX.log"};
        String[] trimmedFileTemplates = new String[]{"{geno}Geno_{pheno}_maf_" + MAF_THRESHOLD + ".gz", "{geno}Geno_{pheno}_maf_" + MAF_THRESHOLD + "-chrX.gz"};
        String[] chrLabels = new String[]{"1:22", "X"};

        String trimmedLabel = "maf_" + MAF_THRESHOLD;
        ArrayList<File[]> rowBoltFilesToTrim = new ArrayList<>();

        try (SimpleFileWriter incompleteWriter = new SimpleFileWriter(incompleteFileList, false)) {

            incompleteWriter.writeLine("# Incomplete analyes" + Instant.now());
            incompleteWriter.writeLine("Phenotype", "Individual", "Chromosome");

            try (SimpleFileWriter failedWriter = new SimpleFileWriter(failedFileList, false)) {

                incompleteWriter.writeLine("# Failed analyses" + Instant.now());
                incompleteWriter.writeLine("Phenotype", "Individual", "Chromosome");

                for (File phenoFolder : boltResultsFolder.listFiles()) {

                    if (phenoFolder.isDirectory()) {

                        String pheno = phenoFolder.getName();

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
                                String trimmedFileName = trimmedFileTemplates[chrI]
                                        .replace("{geno}", geno)
                                        .replace("{pheno}", pheno);

                                File logFile = new File(phenoFolder, logFileName);
                                File boltFile = new File(phenoFolder, boltFileName);
                                File boltTabFile = new File(phenoFolder, boltTabFileName);

                                if (logFile.exists()) {

                                    if (boltFile.exists() && boltTabFile.exists()) {

                                        File trimmedFolder = new File(phenoFolder, trimmedLabel);

                                        if (!trimmedFolder.exists()) {

                                            trimmedFolder.mkdir();

                                        }

                                        File trimmedFile = new File(trimmedFolder, trimmedFileName);

                                        if (true || !trimmedFile.exists()) {

                                            rowBoltFilesToTrim.add(
                                                    new File[]{boltFile, trimmedFile}
                                            );

                                            if (report.length() > 0) {

                                                report.append(", ");

                                            }
                                            report.append(geno + " chr" + chrLabel + " to trim");

                                        } else {

                                            if (report.length() > 0) {

                                                report.append(", ");

                                            }
                                            report.append(geno + " chr" + chrLabel + " done");

                                        }

                                    } else {

                                        failedWriter.writeLine(pheno, geno, chrLabel);

                                        if (report.length() > 0) {

                                            report.append(", ");

                                        }
                                        report.append(geno + " chr" + chrLabel + " failed");

                                        logFile.delete();
                                        boltFile.delete();
                                        boltTabFile.delete();

                                    }

                                } else {

                                    incompleteWriter.writeLine(pheno, geno, chrLabel);

                                    if (report.length() > 0) {

                                        report.append(", ");

                                    }
                                    report.append(geno + " chr" + chrLabel + " missing");

                                }
                            }
                        }

                        Instant end = Instant.now();
                        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                        System.out.println(Instant.now() + "    " + pheno + ": " + report + " (" + durationSeconds + " s)");

                    }
                }
            }
        }

        // Trim bolt results by maf and info score
        rowBoltFilesToTrim.stream()
                .parallel()
                .forEach(
                        filePair -> processBoltFile(
                                filePair[0],
                                filePair[1],
                                rowBoltFilesToTrim.size()
                        )
                );

    }

    private static void processBoltFile(
            File boltFile,
            File trimmedFile,
            int nFiles
    ) {

        System.out.println(Instant.now() + "    Trimming " + boltFile + " (" + trimmingProgress++ + " of " + nFiles + ").");

        Instant begin = Instant.now();

        int nVariants = 0;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(boltFile, false)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(trimmedFile, true)) {

                String line = reader.readLine();
                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    if (!lineSplit[11].equals("-nan") && !lineSplit[11].equals("nan")) {

                        double maf = Double.parseDouble(lineSplit[6]);

                        if (maf >= MAF_THRESHOLD && maf <= 1.0 - MAF_THRESHOLD) {

                            boolean summaryStatsOK = true;

                            try {

                                double beta = Double.parseDouble(lineSplit[10]);

                                if (Double.isNaN(beta) || Double.isInfinite(beta)) {

                                    summaryStatsOK = false;

                                }

                                double se = Double.parseDouble(lineSplit[11]);

                                if (Double.isNaN(se) || Double.isInfinite(se)) {

                                    summaryStatsOK = false;

                                }

                                double p = Double.parseDouble(lineSplit[15]);

                                if (Double.isNaN(p) || Double.isInfinite(p)) {

                                    summaryStatsOK = false;

                                }

                            } catch (Exception e) {

                                summaryStatsOK = false;

                            }

                            if (summaryStatsOK) {

                                writer.writeLine(line);

                                nVariants++;

                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            
            trimmedFile.delete();
            
        }

        Instant end = Instant.now();
        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Trimming " + boltFile + " done, " + nVariants + " variants remaining (" + durationSeconds + " s)");

    }
}
