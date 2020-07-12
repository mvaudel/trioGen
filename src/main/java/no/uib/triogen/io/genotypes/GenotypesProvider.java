package no.uib.triogen.io.genotypes;

import no.uib.triogen.model.family.ChildToParentMap;

/**
 * The GenotypesProvider provides the genotypes for a given variant.
 *
 * @author Marc Vaudel
 */
public interface GenotypesProvider {
    
    /**
     * Parses the genotypes.
     * 
     * @param childToParentMap The child to parent map to use.
     */
    public void parse(
            ChildToParentMap childToParentMap
    );
    
    /**
     * Returns the hard-called genotype for a given sample. 0: 0|0 1: 1|0 2: 0|1 3: 1|1
     *
     * @param sampleId The id of the sample.
     *
     * @return The hard-called genotype.
     */
    public short getGenotype(
            String sampleId
    );
    
    /**
     * Returns the number of alternative alleles using hard calls.
     *
     * @param sampleId The id of the sample.
     *
     * @return The hard-called number of alternative alleles.
     */
    public short getNAlt(
            String sampleId
    );
    
    /**
     * Returns the dosages for the given sample as an array [p0, p1, p2].
     *
     * @param sampleId The id of the sample.
     *
     * @return The dosages.
     */
    public float[] getDosages(
            String sampleId
    );
    
    /**
     * Returns the number of alternative alleles using dosages.
     *
     * @param sampleId The id of the sample.
     *
     * @return The number of alternative alleles using dosages.
     */
    public double getNAltDosages(
            String sampleId
    );
    
    /**
     * Returns the number of alternative alleles for each h of a trio according to the nomenclature of Chen et al. (https://doi.org/10.1101/737106).
     *
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     *
     * @return An array containing the number of alternative alleles.
     */
    public short[] getNAltH(
            String childId,
            String motherId,
            String fatherId
    );

    /**
     * Returns the id of the variant.
     
     * @return the id of the variant
     */
    public String getVariantID();

    /**
     * Returns the contig of the variant.
     
     * @return the contig of the variant
     */
    public String getContig();

    /**
     * Returns the bp of the variant.
     
     * @return the bp of the variant
     */
    public int getBp();

    /**
     * Returns the reference allele of the variant.
     
     * @return the reference allele of the variant
     */
    public String getRef();

    /**
     * Returns the alternative allele of the variant.
     *
     * @return the alternative allele of the variant
     */
    public String getAlt();
    
    /**
     * Returns a boolean indicating whether the marker is genotyped.
     *
     * @return A boolean indicating whether the marker is genotyped.
     */
    public boolean genotyped();
    
    /**
     * Sets the p0s for all parents of the given child ids.
     * 
     * @param childIds The ids of the children.
     * @param childToParentMap The child to parent map.
     */
    public void setParentP0s(
            String[] childIds,
            ChildToParentMap childToParentMap
    );
    
    /**
     * Returns the p0s for all parents of the given child ids. Mothers first and then fathers, in the order of the child ids.
     * 
     * @return The p0s in an array.
     */
    public float[] getParentsDosageP0sCache();
    
    /**
     * Returns the sum of p0s for all parents of the given child ids.
     * 
     * @return The sum of p0s.
     */
    public double getParentsDosageP0Cache();
    
    /**
     * Returns the p0s for all parents of the given child ids. Mothers first and then fathers, in the order of the child ids.
     * 
     * @return The p0s in an array.
     */
    public boolean[] getParentsGenotypeP0sCache();
    
    /**
     * Returns the sum of p0s for all parents of the given child ids.
     * 
     * @return The sum of p0s.
     */
    public int getParentsGenotypeP0Cache();
    
    /**
     * Empty the caches for genotypes and dosages.
     */
    public void emptyGenotypeDosageCaches();

}
