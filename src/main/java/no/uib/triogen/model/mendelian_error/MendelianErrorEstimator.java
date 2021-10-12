package no.uib.triogen.model.mendelian_error;

import no.uib.triogen.io.genotypes.bgen.variant_data.BgenVariantTrioData;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Convenience class to estimate the prevalence of Mendelian errors.
 *
 * @author Marc Vaudel
 */
public class MendelianErrorEstimator {

    /**
     * Estimates the prevalence of Mendelian errors by comparing the number of
     * 001*, *100, 110*, *011 resulting in +2 or -1 tested alleles for the
     * parents compared to the number of such trios expected according to the
     * maf. The check is only done for diploid children.
     *
     * @param variantData The bgen data on this variant.
     * @param childToParentMap The child to parent map to use.
     * @param testedAlleleIndex The index of the tested allele.
     *
     * @return The prevalence as a ratio between the number of observed errors
     * compared to the number of expected trios.
     */
    public static double estimateMendelianErrorPrevalence(
            BgenVariantTrioData variantData,
            ChildToParentMap childToParentMap,
            int testedAlleleIndex
    ) {

        double minusOne = 0;
        double two = 0;
        int nDiploidChildren = 0;

        for (String childId : childToParentMap.children) {

            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            if (variantData.contains(childId) && variantData.getPloidy(childId) == 2) {

                nDiploidChildren++;

                double[] hs = variantData.getHaplotypes(childId, motherId, fatherId, testedAlleleIndex);

                if (hs[0] <= -0.5) {

                    minusOne += 1;

                }
                if (hs[3] <= -0.5) {

                    minusOne += 1;

                }
                if (hs[0] >= 1.5) {

                    two += 1;

                }
                if (hs[3] >= 1.5) {

                    two += 1;

                }
            }
        }

        if (nDiploidChildren == 0) {

            return 0;

        }

        double alleleFrequency = variantData.getAlleleFrequency(testedAlleleIndex);

        double possibleMinusOne = 2.0 * (1 - alleleFrequency) * (1 - alleleFrequency) * alleleFrequency * nDiploidChildren; // number of trios with 001* *100
        double possibleTwo = 2.0 * alleleFrequency * alleleFrequency * (1 - alleleFrequency) * nDiploidChildren; // number of trios with 110* *011

        double measurableErrors = possibleMinusOne + possibleTwo;

        if (measurableErrors <= 10.0) { // Disable if less than 10 cases

            return Double.NaN;

        }

        return ((double) (minusOne + two)) / measurableErrors;

    }
}
