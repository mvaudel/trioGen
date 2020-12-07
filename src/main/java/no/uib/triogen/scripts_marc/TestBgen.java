package no.uib.triogen.scripts_marc;

import java.io.File;
import no.uib.triogen.io.genotypes.bgen.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.BgenUtils;

/**
 *
 *
 * @author Marc Vaudel
 */
public class TestBgen {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        File bgenFile = new File("C:\\Github\\trioGen\\tmp\\X.bgen");
        
        try {
        
        BgenIndex.getBgenIndex(bgenFile);
        
        } catch (Exception e) {
            
            e.printStackTrace();
            
        }
        
    }

}
