package no.uib.triogen.scripts_marc.bolt;

import no.uib.triogen.scripts_marc.mobarun.*;
import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * This script Cleans the organization of BOLT-LMM files on MoBa.
 *
 * @author Marc Vaudel
 */
public class BoltCleanUp {

    /**
     * Main method.
     *
     * Three arguments expected: fileIn fileOut mafThreshold
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // Get file paths
        File fileIn = new File(args[0]);
        File fileOut = new File(args[1]);
        double mafThreshold = Double.parseDouble(args[2]);

        System.out.println(Instant.now() + "    Filtering " + fileIn + " to " + fileOut + " with maf threshold of " + mafThreshold + ".");

        Instant begin = Instant.now();

        int nVariants = 0;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(fileIn, false)) {

            try (SimpleFileWriter commonWriter = new SimpleFileWriter(fileOut, true)) {

                String line = reader.readLine();
                commonWriter.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    if (!lineSplit[11].equals("-nan") && !lineSplit[11].equals("nan")) {

                        double maf = Double.parseDouble(lineSplit[6]);

                        boolean common = false;

                        if (maf >= mafThreshold && maf <= 1.0 - mafThreshold) {

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

                            }
                        }
                    }
                }
            }
        }

        Instant end = Instant.now();
        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Splitting " + fileIn + " done, " + nVariants + " common variants (" + durationSeconds + " s)");

    }
}
