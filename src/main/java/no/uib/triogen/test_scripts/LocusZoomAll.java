package no.uib.triogen.test_scripts;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.export.LocusZoomExtractor;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts locus zoom data.
 *
 * @author Marc Vaudel
 */
public class LocusZoomAll {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File shFile = new File("scripts/test_moba_lz.sh");

            try ( SimpleFileWriter writer = new SimpleFileWriter(shFile, false)) {

                writer.writeLine("#!/usr/bin/env bash");
                writer.newLine();

                long start = Instant.now().getEpochSecond();

                System.out.println("  " + Instant.now() + " Building locus zoom plots.");

//                int[] targetChr = new int[]{1, 10, 2, 3, 5, 6, 7, 9};
                int[] targetChr = new int[]{3};
                String[] phenos = new String[]{"z_bmi0", "z_bmi1", "z_bmi2", "z_bmi3", "z_bmi4", "z_bmi5", "z_bmi6", "z_bmi7", "z_bmi8", "z_bmi9", "z_bmi10", "z_bmi11"};

                for (int chr : targetChr) {

                    File targetFile = new File("resources/targets_chr" + chr);
                    File resultFile = new File("docs/lm_test/target/" + chr + ".lm_target.gz");

                    try ( SimpleFileReader reader = SimpleFileReader.getFileReader(targetFile)) {

                        String targetLine;
                        while ((targetLine = reader.readLine()) != null) {

                            if (!targetLine.startsWith("#") && !targetLine.startsWith("id")) {

                                String[] lineSplit = targetLine.split("\t");

                                String variantId = lineSplit[0];

                                File ldFile = new File("docs/ld_matrix_test/dos_chr_" + chr + "_" + variantId + ".tld");

                                for (String phenoName : phenos) {

                                    String outputPath = "docs/lz/data/" + variantId + "_" + phenoName + "_locusZoomData.gz";
                                    File outputFile = new File(outputPath);

                                    String genesPath = "docs/lz/data/" + variantId + "_" + phenoName + "_locusZoomGenes.gz";
                                    File genesFile = new File(genesPath);

                                    System.out.println("      " + Instant.now() + " Extracting data for locus zoom of " + variantId + " " + phenoName + ".");

                                    LocusZoomExtractor.writeData(
                                            phenoName,
                                            variantId,
                                            1000000,
                                            37,
                                            resultFile,
                                            ldFile,
                                            outputFile,
                                            genesFile
                                    );

                                    writer.writeLine("Rscript src/R/locusZoom.R " + outputPath + " " + genesPath + " " + variantId + " " + phenoName + " docs/lz/plots/" + variantId + "_" + phenoName + "_locusZoom ~/R");

                                }
                            }
                        }
                    }
                }

                long end = Instant.now().getEpochSecond();
                long duration = end - start;

                System.out.println("  " + Instant.now() + " Building locus zoom plots completed (" + duration + " s).");

            }
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
