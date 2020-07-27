package no.uib.triogen.model.mendelian_error;

import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.maf.MafEstimator;

/**
 * Convenience class to estimate the prevalence of Mendelian errors.
 *
 * @author Marc Vaudel
 */
public class MendelianErrorEstimator {

    /**
     * Estimates the prevalence of Mendelian errors by comparing the number of 001*, *100, 110*, *011 resulting in +2 or -1 alternative alleles for the parents compared to the number of such trios expected according to the maf.
     * 
     * @param genotypesProvider The genotypes provider for the variant to inspect.
     * @param childToParentMap The child to parent map to use.
     * 
     * @return The prevalence as a ratio between the number of observed errors compared to the number of expected trios.
     */
    public static double estimateMendelianErrorPrevalence(
            GenotypesProvider genotypesProvider,
            ChildToParentMap childToParentMap
    ) {

        double maf = MafEstimator.getMaf(genotypesProvider, childToParentMap);

        double minusOne = 0.0;
        double two = 0.0;

        for (String childId : childToParentMap.children) {

            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            short[] hs = genotypesProvider.getNAltH(childId, motherId, fatherId);

            if (hs[1] == -1 || hs[3] == -1) {

                minusOne += 1.0;

            } else if (hs[1] == 2 || hs[3] == 2) {

                two += 1.0;

            }
        }

        double expectedMinusOne = childToParentMap.children.length * 2 * (1 - maf) * (1 - maf) * maf; // number of trios with 001* *100
        double expectedTwo = childToParentMap.children.length * 2 * maf * maf * (1 - maf); // number of trios with 110* *011

        return (minusOne + two) / (expectedMinusOne + expectedTwo);

    }

}
