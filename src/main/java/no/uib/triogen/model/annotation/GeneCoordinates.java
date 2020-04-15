package no.uib.triogen.model.annotation;

/**
 * Gene coordinates as retrieved from the Ensembl API.
 *
 * @author Marc Vaudel
 */
public class GeneCoordinates {

    public final String biotype;
    public final String name;
    public final int start;
    public final int end;

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
