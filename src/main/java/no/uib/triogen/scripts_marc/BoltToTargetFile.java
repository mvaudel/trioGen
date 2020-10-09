package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * This script generates a target file from a Bolt results file
 *
 * @author Marc Vaudel
 */
public class BoltToTargetFile {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File prunedFile = new File(args[0]);

        System.out.println(Instant.now() + "    Processing " + prunedFile.getName());

        Instant begin = Instant.now();

        File sourceFolder = prunedFile.getParentFile();
        File destinationFolder = new File(sourceFolder, "targets");

        String fileName = prunedFile.getName();

        String genoPheno = fileName.substring(0, fileName.length() - 21);

        String chrXFileName = fileName.substring(0, fileName.length() - 3) + "-chrX.gz";
        String targetFileName = fileName.substring(0, fileName.length() - 3) + "_targets";

        File destinationFile = new File(destinationFolder, targetFileName);

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, false)) {

            writer.writeLine("snp", "chr", "start", "end", "source");

            File[] files = new File[]{prunedFile, new File(sourceFolder, chrXFileName)};

            for (File file : files) {

                if (file.exists()) {

                    try (SimpleFileReader reader = SimpleFileReader.getFileReader(prunedFile)) {

                        String line = reader.readLine();

                        while ((line = reader.readLine()) != null) {

                            String[] lineSplit = line.split("\t");

                            double maf = Double.parseDouble(lineSplit[6]);

                            if (maf > 0.5) {

                                maf = 1 - maf;

                            }

                            if (maf >= 0.005) {

                                String snp = lineSplit[0];
                                String chr = lineSplit[1];
                                String bp = lineSplit[2];

                                writer.writeLine(snp, chr, bp, bp, genoPheno);

                            }
                        }
                    }
                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Processing " + prunedFile.getName() + " finished (" + durationSeconds + " s)");

    }
}
