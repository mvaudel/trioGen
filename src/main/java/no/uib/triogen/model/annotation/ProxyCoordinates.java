package no.uib.triogen.model.annotation;

/**
 * Proxy coordinates as returned from Ensembl.
 *
 * @author Marc Vaudel
 */
public class ProxyCoordinates {
    
    /**
     * The identifier.
     */
    public final String proxySnp;
    /**
     * The contig.
     */
    public final String contig;
    /**
     * The start position.
     */
    public final int start;
    /**
     * The end position.
     */
    public final int end;
    /**
     * r2.
     */
    public final double r2;
    
    /**
     * Constructor.
     * 
     * @param proxySnp The identifier of the proxy.
     * @param contig The contig.
     * @param start The start position of the variant.
     * @param end The end position of the variant.
     * @param r2 The r2 with the queried variant.
     */
    public ProxyCoordinates(
            String proxySnp,
            String contig,
            int start,
            int end,
            double r2
    ) {
        
        this.proxySnp = proxySnp;
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.r2 = r2;
        
    }

}
