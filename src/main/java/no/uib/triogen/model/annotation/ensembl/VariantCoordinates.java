package no.uib.triogen.model.annotation.ensembl;

/**
 * Variant coordinates as returned from Ensembl.
 *
 * @author Marc Vaudel
 */
public class VariantCoordinates {
    
    /**
     * The contig.
     */
    public final String contig;
    /**
     * The position.
     */
    public final int position;
    
    /**
     * Constructor.
     * 
     * @param contig The contig.
     * @param position The position.
     */
    public VariantCoordinates(
            String contig,
            int position
    ) {
        
        this.contig = contig;
        this.position = position;
        
    }

}
