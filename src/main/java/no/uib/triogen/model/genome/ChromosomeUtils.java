package no.uib.triogen.model.genome;

import java.util.HashMap;

/**
 * Convenience data and functions for chromosomes.
 *
 * @author Marc Vaudel
 */
public class ChromosomeUtils {
    
    /**
     * Map of the chromosome length in BP in GRCh37.p13 (hg19) from Ensembl.
     */
    public static final HashMap<String, Integer> chromosomeLength37 = getchromosomeLength37();
    
    /**
     * Returns the map of the chromosome length in BP in GRCh37.p13 (hg19) from Ensembl.
     * 
     * @return The map of the chromosome length in BP in GRCh37.p13 (hg19) from Ensembl.
     */
    private static HashMap<String, Integer> getchromosomeLength37() {
        
        HashMap<String, Integer> result = new HashMap<>(23);
        
        String[] chrNames = new String[]{"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"};
        int[] chrLength = new int[]{249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560};
        
        for (int i = 0 ; i < chrNames.length ; i++) {
            
            result.put(chrNames[i], chrLength[i]);
            
        }
        
        return result;
        
    }

}
