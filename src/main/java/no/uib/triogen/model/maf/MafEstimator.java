package no.uib.triogen.model.maf;

import no.uib.triogen.io.genotypes.vcf.reader.VcfVariant;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Estimator for the minor allele frequency.
 *
 * @author Marc Vaudel
 */
public class MafEstimator {
    
    /**
     * Returns the maf for the given genotypes provider using all parents of the children in the given child to parent map.
     * 
     * @param genotypesProvider The genotypes provider.
     * @param childToParentMap The child to parent map.
     * 
     * @return the maf.
     */
    public static double getMaf(
            VcfVariant genotypesProvider,
            ChildToParentMap childToParentMap
    ) {
        
        double nAlt = 0.0;
        
        for (String childId : childToParentMap.children) {
            
            String motherId = childToParentMap.getMother(childId);
            float[] dosages = genotypesProvider.getGenotypingProbabilities(motherId);
            nAlt += dosages[1] + 2 * dosages[2];
            
            String fatherId = childToParentMap.getMother(childId);
            dosages = genotypesProvider.getGenotypingProbabilities(fatherId);
            nAlt += dosages[1] + 2 * dosages[2];
            
        }
        
        double nAlleles = (double) (4 * childToParentMap.children.length);
        
        double maf = nAlt / nAlleles;
        
        if (maf > 0.5) {
            
            maf = 1.0 - maf;
            
        }
        
        return maf;
        
    }

}
