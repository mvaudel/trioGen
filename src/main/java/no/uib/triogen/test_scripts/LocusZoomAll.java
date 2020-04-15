package no.uib.triogen.test_scripts;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import no.uib.triogen.export.LocusZoomExtractor;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.ld.LdMatrixReader;

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

            long start = Instant.now().getEpochSecond();

            System.out.println("  " + Instant.now() + " Building locus zoom plots.");

            int[] targetChr = new int[]{1, 10, 2, 3, 5, 6, 7, 9};
            String[] phenos = new String[]{"breastmilk_duration", "formula_freq_6m", "pregnancy_duration", "z_umbilical_chord_length", "z_placenta_weight", "z_bmi0", "z_bmi1", "z_bmi2", "z_bmi3", "z_bmi4", "z_bmi5", "z_bmi6", "z_bmi7", "z_bmi8", "z_bmi9", "z_bmi10", "z_bmi11", "z_mother_height", "z_father_bmi"};

            for (int chr : targetChr) {

                File targetFile = new File("resources/targets_chr" + chr);
                File resultFile = new File("/mnt/work/marc/moba/test_TrioGen/results/chr_" + chr + ".gz");
                File ldFile = new File("/mnt/work/marc/moba/test_TrioGen/ld/chr_" + chr + ".tld");

                try ( SimpleFileReader reader = SimpleFileReader.getFileReader(targetFile)) {

                    String targetLine;
                    while ((targetLine = reader.readLine()) != null) {

                        if (!targetLine.startsWith("#") && !targetLine.startsWith("id")) {

                            String[] lineSplit = targetLine.split("\t");

                            String variantId = lineSplit[0];

                            for (String phenoName : phenos) {

                                File outputFile = new File("docs/lz/" + variantId + "_" + phenoName + "_locusZoomData.gz");
                                File genesFile = new File("docs/lz/" + variantId + "_" + phenoName + "_locusZoomGenes.gz");

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

                                ArrayList<String> rCommand = new ArrayList<>();
                                rCommand.add("Rscript");
                                rCommand.add(outputFile.getAbsolutePath());
                                rCommand.add(genesFile.getAbsolutePath());
                                rCommand.add(variantId);
                                rCommand.add(phenoName);
                                rCommand.add("docs/lz/" + variantId + "_" + phenoName);
                                rCommand.add("~/R");

                                System.out.println("      " + Instant.now() + " Running LousZoom R script.");
                                System.out.println("        " + rCommand.stream().collect(Collectors.joining(" ")));

                                ProcessBuilder pb = new ProcessBuilder(rCommand);
                                Process p = pb.start();
                                p.waitFor();

                            }
                        }
                    }
                }
            }

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            System.out.println("  " + Instant.now() + " Building locus zoom plots completed (" + duration + " s).");

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
