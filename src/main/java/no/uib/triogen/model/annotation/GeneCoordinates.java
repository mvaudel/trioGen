package no.uib.triogen.model.annotation;

/**
 * Gene coordinates as retrieved from the Ensembl API.
 *
 * @author Marc Vaudel
 */
public class GeneCoordinates {

    /**
     * The biotype.
     */
    public final String biotype;
    /**
     * The name.
     */
    public final String name;
    /**
     * The start position.
     */
    public final int start;
    /**
     * The end position.
     */
    public final int end;

    /**
     * Constructor.
     * 
     * @param biotype The biotype.
     * @param name The name.
     * @param start The start position.
     * @param end The end position.
     */
    public GeneCoordinates(
            String biotype,
            String name,
            int start,
            int end
    ) {
        this.biotype = biotype;
        this.name = name;
        this.start = start;
        this.end = end;
    }

}
