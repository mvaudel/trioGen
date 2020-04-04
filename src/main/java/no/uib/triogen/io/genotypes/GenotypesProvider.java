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
     */
    public void parse();
    
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
     * Returns the hs for a trio according to the nomenclature of Chen et al. (https://doi.org/10.1101/737106).
     *
     * @param sampleId The id of the child in the trio.
     * @param childToParentMap The child to parent map.
     *
     * @return An array containing the hs.
     */
    public short[] getH(
            ChildToParentMap childToParentMap,
            String sampleId
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

}
