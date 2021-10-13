package no.uib.triogen.scripts_marc.phasing;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Reformats the hapmap_recombination_2011-01_phaseII_B37 genetic maps for
 * shapeit.
 *
 * @author Marc Vaudel
 */
public class ReformatGeneticMaps {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String hapmapGenericPath = "/home/marcvaue/resources/hapmap_recombination_2011-01_phaseII_B37/genetic_map_GRCh37_chr{chr}.txt"; // QC node
        String outputGenericPath = "/home/marcvaue/resources/hapmap_recombination_2011-01_phaseII_B37/genetic_map_GRCh37_chr{chr}_shapeit.gz"; // QC node

        for (String chromosme : getChromosomes()) {

            System.out.println(Instant.now() + " - Processing chromosome " + chromosme);

            String hapmapFile = hapmapGenericPath.replace("{chr}", chromosme);
            String outputFile = outputGenericPath.replace("{chr}", chromosme);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(new File(hapmapFile))) {

                try (SimpleFileWriter writer = new SimpleFileWriter(new File(outputFile), true)) {

                    writer.writeLine("pposition rrate gposition");

                    String line = reader.readLine();

                    while ((line = reader.readLine()) != null) {

                        String[] lineSplit = line.split("\t");

                        String newLine = String.join(" ", lineSplit[1], lineSplit[2], lineSplit[3]);

                        writer.writeLine(newLine);

                    }
                }
            }
        }
    }

    private static String[] getChromosomes() {

        String[] chromosomes = new String[25];

        for (int chr = 1; chr <= 22; chr++) {

            chromosomes[chr - 1] = Integer.toString(chr);

        }

        chromosomes[22] = "X";
        chromosomes[23] = "chrX_par1";
        chromosomes[24] = "chrX_par2";

        return chromosomes;

    }

}
