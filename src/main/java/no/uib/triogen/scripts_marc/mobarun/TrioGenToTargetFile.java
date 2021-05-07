package no.uib.triogen.scripts_marc.mobarun;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;

/**
 * This script generates a target file from a Bolt results file
 *
 * @author Marc Vaudel
 */
public class TrioGenToTargetFile {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File resultsFolder = new File(args[0]);

        String pheno = resultsFolder.getName();

        System.out.println(Instant.now() + "    Processing " + pheno);

        Instant begin = Instant.now();

        File prunedFolder = new File(resultsFolder, "pruned");
        File destinationFolder = new File(prunedFolder, "targets");

        String targetFileName = pheno + "_targets";

        File destinationFile = new File(destinationFolder, targetFileName);

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, false)) {

            writer.writeLine("snp", "chr", "start", "end", "source");

            for (File file : prunedFolder.listFiles()) {

                if (!file.isDirectory() && file.getName().endsWith("_pruned.gz")) {

                    String[] nameSplit = file.getName().split("\\.");

                    StringBuilder sb = new StringBuilder();

                    for (int i = 2; i < nameSplit.length - 2; i++) {

                        if (i > 2) {

                            sb.append('.');

                        }

                        sb.append(nameSplit[i]);

                    }

                    String pName = sb.toString();

                    String source = String.join("_", pheno, pName);

                    try (SimpleFileReader reader = SimpleFileReader.getFileReader(file)) {

                        String line = reader.readLine();

                        while ((line = reader.readLine()) != null) {

                            String[] lineSplit = line.split("\t");

                            String snp = lineSplit[3];
                            String chr = lineSplit[1];
                            String bp = lineSplit[2];

                            writer.writeLine(snp, chr, bp, bp, source);

                        }
                    }
                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Processing " + pheno + " finished (" + durationSeconds + " s)");

    }
}
