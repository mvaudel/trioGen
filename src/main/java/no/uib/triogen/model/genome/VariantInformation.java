package no.uib.triogen.model.genome;

/**
 * 
 *
 * @author Marc Vaudel
 */
public class VariantInformation {
    
    /**
     * The identifier of the variant.
     */
    public final String id;
    
    /**
     * The rsId of the variant.
     */
    public final String rsId;
    
    /**
     * The contig of the variant.
     */
    public final String contig;
    
    /**
     * The location of the variant.
     */
    public final int bp;
    
    /**
     * The alleles.
     */
    public final String[] alleles;
    
    /**
     * Constructor.
     * 
     * @param id The identifier of the variant.
     * @param rsId The rsid of the variant.
     * @param contig The contig of the variant.
     * @param bp The location of the variant.
     * @param alleles The alleles.
     */
    public VariantInformation(
            String id,
            String rsId,
            String contig,
            int bp,
            String[] alleles
    ) {
        
        this.id = id;
        this.rsId = rsId;
        this.contig = contig;
        this.bp = bp;
        this.alleles = alleles;
        
    }
}
