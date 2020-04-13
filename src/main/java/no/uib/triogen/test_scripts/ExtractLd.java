package no.uib.triogen.test_scripts;

import java.io.File;
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

            File ldFile = new File("docs/tmp/chr_9.tld");
            LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);

            String[] variantIds = ldMatrixReader.variantIds;

            for (String variantId : variantIds) {

                HashMap<String, Double> variantLdMap = ldMatrixReader.getR2(variantId);

                if (variantLdMap != null) {

                    for (Entry<String, Double> entry : variantLdMap.entrySet()) {

                        if (entry.getValue() < 0.0 || entry.getValue() > 1.0) {
                            
                            int debug = 1;
                        
                        }
                    }
                }
            }

            int debug = 1;

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
