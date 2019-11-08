package no.uib.triogen.io.genotypes;

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
     * Returns the genotype for a given sample. 0: 0|0 1: 1|0 2: 0|1 3: 1|1
     *
     * @param sampleId the id of the sample
     *
     * @return the genotype
     */
    public int getGenotype(
            String sampleId
    );

    /**
     * Returns the id of the variant.
     
     * @return the id of the variant
     */
    public String getVariantID();

}
