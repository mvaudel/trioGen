package no.uib.triogen.scripts_marc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.model.ld.R2;

/**
 * Tests the compression and decompression methods.
 *
 * @author Marc Vaudel
 */
public class CompressionTest {

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

                ArrayList<R2> variantLdMap = ldMatrixReader.getR2(variantId);

                if (variantLdMap != null) {

                    for (R2 r2 : variantLdMap) {

                        if (r2.r2Value < 0.0 || r2.r2Value > 1.0) {
                            
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
