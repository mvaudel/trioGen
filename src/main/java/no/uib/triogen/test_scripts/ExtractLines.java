package no.uib.triogen.test_scripts;

import java.io.File;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts specific lines of a file.
 *
 * @author Marc Vaudel
 */
public class ExtractLines {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        int start = 0;
        int end = 100;

        try {

            File gzFile = new File("C:\\Projects\\TrioGen\\vcf\\ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.2.vcf.gz");
            File outFile = new File("C:\\Projects\\TrioGen\\vcf\\test.txt");
    
            try (SimpleFileReader reader = SimpleFileReader.getFileReader(gzFile)) {
                
                try (SimpleFileWriter writer = new SimpleFileWriter(outFile, false)) {
                    
                    int count = 0;
                    
                    String line;
                    while((line = reader.readLine()) != null && ++count <= end) {
                        
                        if (count >= start) {
                            
                            writer.writeLine(line);
                            
                        }
                    }
                }
            }
    
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
