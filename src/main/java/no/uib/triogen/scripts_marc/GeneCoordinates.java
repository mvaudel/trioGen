package no.uib.triogen.scripts_marc;

import java.io.File;
import java.util.ArrayList;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.model.annotation.EnsemblAPI;

/**
 *
 *
 * @author Marc Vaudel
 */
public class GeneCoordinates {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        int buildNumber = 37;
        String targetContig = "11";
        int bpStart = 1629111;
        int bpEnd = 3349109;
        
        File geneFile = new File("C:\\Github\\placenta_weight\\src\\IGF2_KCNQ1\\resources\\genes.gz");
        
        try {

                try ( SimpleFileWriter writer = new SimpleFileWriter(geneFile, true)) {

                    String ensemblVersion = EnsemblAPI.getEnsemblVersion(buildNumber);

                    writer.writeLine("# Ensembl version: " + ensemblVersion);
                    writer.writeLine("biotype", "name", "start", "end");

                    ArrayList<no.uib.triogen.model.annotation.GeneCoordinates> geneCoordinatesList = EnsemblAPI.getGeneCoordinates(
                            targetContig,
                            bpStart,
                            bpEnd,
                            buildNumber
                    );

                    for (no.uib.triogen.model.annotation.GeneCoordinates geneCoordinates : geneCoordinatesList) {

                        writer.writeLine(
                                geneCoordinates.biotype,
                                geneCoordinates.name,
                                Integer.toString(geneCoordinates.start),
                                Integer.toString(geneCoordinates.end)
                        );
                    }
                }
        
        } catch (Exception e) {
            
            e.printStackTrace();
            
        }
        
    }
}
