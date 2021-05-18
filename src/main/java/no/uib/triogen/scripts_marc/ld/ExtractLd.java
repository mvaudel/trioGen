package no.uib.triogen.scripts_marc.ld;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.model.ld.R2;

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

            File ldFile = new File("tmp/20.tld");
            LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);
            
            
                ArrayList<R2> result = ldMatrixReader.getR2("rs6040450");

            String[] variantIds = ldMatrixReader.variantIds;

            for (String variantId : variantIds) {
                
                if (variantId.equals("20_11209782_C_T")) {
                    int debug = 1;
                }

                ArrayList<R2> variantLdMap = ldMatrixReader.getR2(variantId);

                if (variantLdMap != null) {

                    for (R2 r2 : variantLdMap) {

                        if (r2.r2Value < -0.001 || r2.r2Value > 1.001) {

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
