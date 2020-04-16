package no.uib.triogen.test_scripts;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import no.uib.triogen.io.ld.LdMatrixReader;

/**
 * Extracts specific lines of a file.
 *
 * @author Marc Vaudel
 */
public class ExtractLd {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File ldFile = new File("docs/tmp/ch_chr_7.tld");
            LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);

            String[] variantIds = ldMatrixReader.variantIds;

            for (String variantId : variantIds) {

                HashMap<String, Double> variantLdMap = ldMatrixReader.getR2(variantId);

                if (variantLdMap != null) {

                    for (Entry<String, Double> entry : variantLdMap.entrySet()) {

                        if (entry.getValue() < -0.001 || entry.getValue() > 1.001) {

                            throw new IllegalArgumentException("Incorrect r2");

                        }
                    }
                }
            }

            Instant begin = Instant.now();

            Arrays.stream(variantIds)
                    .forEach(
                            variantId -> ldMatrixReader.getR2(variantId)
                    );

            Instant end = Instant.now();

            long timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("LD matrix iterated in " + timeInSec + "s using a single thread.");

            begin = Instant.now();

            Arrays.stream(variantIds)
                    .parallel()
                    .forEach(
                            variantId -> ldMatrixReader.getR2(variantId)
                    );

            end = Instant.now();

            timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("LD matrix iterated in " + timeInSec + "s using multiple threads.");

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
