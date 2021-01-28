package no.uib.triogen.model.mendelian_error;

import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Convenience class to estimate the prevalence of Mendelian errors.
 *
 * @author Marc Vaudel
 */
public class MendelianErrorEstimator {

    /**
     * Estimates the prevalence of Mendelian errors by comparing the number of 001*, *100, 110*, *011 resulting in +2 or -1 tested alleles for the parents compared to the number of such trios expected according to the maf.
     * 
     * @param genotypesProvider The genotypes provider for the variant to inspect.
     * @param childToParentMap The child to parent map to use.
     * @param testedAlleleIndex The index of the tested allele.
     * 
     * @return The prevalence as a ratio between the number of observed errors compared to the number of expected trios.
     */
    public static double estimateMendelianErrorPrevalence(
            BgenVariantData genotypesProvider,
            ChildToParentMap childToParentMap,
            int testedAlleleIndex
    ) {

        double minusOne = 0;
        double two = 0;

        for (String childId : childToParentMap.children) {

            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            double[] hs = genotypesProvider.getHaplotypes(childId, motherId, fatherId, testedAlleleIndex);

            if (hs[1] <= -0.5 || hs[3] <= -0.5) {

                minusOne += 1;

            } 
            if (hs[1] >= 1.5 || hs[3] >= 1.5) {

                two += 1;

            }
        }

        double alleleFrequency = genotypesProvider.getAlleleFrequency(testedAlleleIndex);

        double expectedMinusOne = childToParentMap.children.length * 2 * (1 - alleleFrequency) * (1 - alleleFrequency) * alleleFrequency; // number of trios with 001* *100
        double expectedTwo = childToParentMap.children.length * 2 * alleleFrequency * alleleFrequency * (1 - alleleFrequency); // number of trios with 110* *011

        return ((double) (minusOne + two)) / (expectedMinusOne + expectedTwo);

    }

}
