package no.uib.triogen.scripts_marc.ld;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.processing.ld.LdMatrixComputer;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Extracts specific lines of a file.
 *
 * @author Marc Vaudel
 */
public class LdTest {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File ldFile = new File("tmp/12_200k.tld");
            
            LdMatrixReader reader = new LdMatrixReader(ldFile);
            
            ArrayList<R2> r2s = reader.getR2("12_66343400_G_C");
            String rsidA = reader.getRsId("12_66343400_G_C");
            
            for (R2 r2 : r2s) {
                
                String id = reader.getId(r2.variantB);
                String rsid = reader.getRsId(id);
                
                if (id.equals("12_66327632_A_G")) {
                    
                    int debug = 1;
                    
                }
                
            }

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
