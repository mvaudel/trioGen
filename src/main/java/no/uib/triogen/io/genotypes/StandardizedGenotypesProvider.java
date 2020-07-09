package no.uib.triogen.io.genotypes;

import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Convenience class to standardize the genotypes.
 *
 * @author Marc Vaudel
 */
public class StandardizedGenotypesProvider {

    /**
     * The genotypes provider to use.
     */
    private final GenotypesProvider genotypesProvider;
    /**
     * Cache for the mean values. 0: h1 1: h2 2: h3 3: h4 4: child number of
     * alternative alleles based on hard calls 5: mother number of alternative
     * alleles based on hard calls 6: father number of alternative alleles based
     * on hard calls 7: child number of alternative alleles based on dosages 8:
     * mother number of alternative alleles based on dosages 9: father number of
     * alternative alleles based on dosages
     */
    private double[] meanCache = new double[10];
    /**
     * Cache for the sd values. 0: h1 1: h2 2: h3 3: h4 4: child number of
     * alternative alleles based on hard calls 5: mother number of alternative
     * alleles based on hard calls 6: father number of alternative alleles based
     * on hard calls 7: child number of alternative alleles based on dosages 8:
     * mother number of alternative alleles based on dosages 9: father number of
     * alternative alleles based on dosages
     */
    private double[] sdCache = new double[10];

    /**
     * Constructor.
     *
     * @param genotypesProvider The genotype provider to use for
     * standardization.
     * @param childToParentMap The child to parent map.
     */
    public StandardizedGenotypesProvider(
            GenotypesProvider genotypesProvider,
            ChildToParentMap childToParentMap
    ) {

        this.genotypesProvider = genotypesProvider;

        int nTrios = childToParentMap.children.length;

        for (String childId : childToParentMap.children) {

            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            short[] nAlth = genotypesProvider.getNAltH(childId, motherId, fatherId);
            short nAltChild = genotypesProvider.getNAlt(childId);
            short nAltMother = genotypesProvider.getNAlt(motherId);
            short nAltFather = genotypesProvider.getNAlt(fatherId);
            double nAltChildDosages = genotypesProvider.getNAltDosages(childId);
            double nAltMotherDosages = genotypesProvider.getNAltDosages(motherId);
            double nAltFatherDosages = genotypesProvider.getNAltDosages(fatherId);

            meanCache[0] += nAlth[0];
            meanCache[1] += nAlth[1];
            meanCache[2] += nAlth[2];
            meanCache[3] += nAlth[3];
            meanCache[4] += nAltChild;
            meanCache[5] += nAltMother;
            meanCache[6] += nAltFather;
            meanCache[7] += nAltChildDosages;
            meanCache[8] += nAltMotherDosages;
            meanCache[9] += nAltFatherDosages;

        }

        for (int i = 0; i < meanCache.length; i++) {

            meanCache[i] = meanCache[i] / nTrios;

        }

        for (String childId : childToParentMap.children) {

            String motherId = childToParentMap.getMother(childId);
            String fatherId = childToParentMap.getFather(childId);

            short[] nAlth = genotypesProvider.getNAltH(childId, motherId, fatherId);
            short nAltChild = genotypesProvider.getNAlt(childId);
            short nAltMother = genotypesProvider.getNAlt(motherId);
            short nAltFather = genotypesProvider.getNAlt(fatherId);
            double nAltChildDosages = genotypesProvider.getNAltDosages(childId);
            double nAltMotherDosages = genotypesProvider.getNAltDosages(motherId);
            double nAltFatherDosages = genotypesProvider.getNAltDosages(fatherId);

            sdCache[0] += (nAlth[0] - meanCache[0]) * (nAlth[0] - meanCache[0]);
            sdCache[1] += (nAlth[1] - meanCache[1]) * (nAlth[1] - meanCache[1]);
            sdCache[2] += (nAlth[2] - meanCache[2]) * (nAlth[2] - meanCache[2]);
            sdCache[3] += (nAlth[3] - meanCache[3]) * (nAlth[3] - meanCache[3]);
            sdCache[4] += (nAltChild - meanCache[4]) * (nAltChild - meanCache[4]);
            sdCache[5] += (nAltMother - meanCache[5]) * (nAltMother - meanCache[5]);
            sdCache[6] += (nAltFather - meanCache[6]) * (nAltFather - meanCache[6]);
            sdCache[7] += (nAltChildDosages - meanCache[7]) * (nAltChildDosages - meanCache[7]);
            sdCache[8] += (nAltMotherDosages - meanCache[8]) * (nAltMotherDosages - meanCache[8]);
            sdCache[9] += (nAltFatherDosages - meanCache[9]) * (nAltFatherDosages - meanCache[9]);

        }

        for (int i = 0; i < meanCache.length; i++) {

            sdCache[i] = Math.sqrt(sdCache[i] / nTrios);

        }
    }

    /**
     * Returns the number of alternative alleles for each h centered and scaled.
     *
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     *
     * @return An array containing the number of alternative alleles centered
     * and scaled.
     */
    public double[] getNAltH(
            String childId,
            String motherId,
            String fatherId
    ) {

        short[] nAlth = genotypesProvider.getNAltH(childId, motherId, fatherId);

        double[] result = new double[4];

        result[0] = (nAlth[0] - meanCache[0]) / sdCache[0];
        result[1] = (nAlth[1] - meanCache[1]) / sdCache[1];
        result[2] = (nAlth[2] - meanCache[2]) / sdCache[2];
        result[3] = (nAlth[3] - meanCache[3]) / sdCache[3];

        return result;

    }

    /**
     * Returns the number of alternative alleles for a child centered and
     * scaled.
     *
     * @param childId The id of the child.
     * @param useDosages If true, dosages are used instead of hard calls.
     *
     * @return The number of alternative alleles for a child centered and
     * scaled.
     */
    public double getNAltChild(
            String childId,
            boolean useDosages
    ) {

        return useDosages
                ? (genotypesProvider.getNAltDosages(childId) - meanCache[7]) / sdCache[7]
                : (genotypesProvider.getNAlt(childId) - meanCache[4]) / sdCache[4];

    }

    /**
     * Returns the number of alternative alleles for a mother centered and
     * scaled.
     *
     * @param motherId The id of the mother.
     * @param useDosages If true, dosages are used instead of hard calls.
     *
     * @return The number of alternative alleles for a mother centered and
     * scaled.
     */
    public double getNAltMother(
            String motherId,
            boolean useDosages
    ) {

        return useDosages
                ? (genotypesProvider.getNAltDosages(motherId) - meanCache[8]) / sdCache[8]
                : (genotypesProvider.getNAlt(motherId) - meanCache[5]) / sdCache[5];

    }

    /**
     * Returns the number of alternative alleles for a father centered and
     * scaled.
     *
     * @param fatherId The id of the father.
     * @param useDosages If true, dosages are used instead of hard calls.
     *
     * @return The number of alternative alleles for a mother centered and
     * scaled.
     */
    public double getNAltFather(
            String fatherId,
            boolean useDosages
    ) {

        return useDosages
                ? (genotypesProvider.getNAltDosages(fatherId) - meanCache[9]) / sdCache[9]
                : (genotypesProvider.getNAlt(fatherId) - meanCache[6]) / sdCache[6];

    }

}
