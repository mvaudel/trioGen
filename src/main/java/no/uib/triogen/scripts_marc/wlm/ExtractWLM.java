package no.uib.triogen.scripts_marc.wlm;

import java.io.File;
import java.time.Instant;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class runs WLM on meta results. Based on work by RN Beaumont.
 *
 * @author Marc Vaudel
 */
public class ExtractWLM {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File wlmFile = new File("C:\\Projects\\placenta_weight\\meta_results\\prs_wlm_200000.gz");
        File destinationFile = new File("C:\\Projects\\placenta_weight\\meta_results\\prs_wlm_200000_pw.gz");
        
        File targetsFile = new File("C:\\Github\\placenta_weight\\resources\\targets\\targets_pw");
        
        int distance = 500000;
        
          VariantList  variantList = VariantList.getVariantList(
                    targetsFile
            );
          
          variantList.index(distance);

        // Load child
        System.out.println(Instant.now() + "    Loading child results.");

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(wlmFile)) {
            
            try(SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {
                
                String line = reader.readLine();
                
                writer.writeLine(line);
                
                while ((line = reader.readLine()) != null) {
                    
                    String[] lineSplit = line.split("\t");
                    
                    String chromosome = lineSplit[2];
                    int bp = Integer.parseInt(lineSplit[3]);
                    
                    if (variantList.include(chromosome, bp)) {
                        
                        writer.writeLine(line);
                        
                    }
                }
            }
        }
    }
}
