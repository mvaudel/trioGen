package no.uib.triogen.io.flat.indexed.gz;

/**
 * Placeholder for the coordinates of a section in the gzip file.
 *
 * @author Marc Vaudel
 */
public class IndexedGzCoordinates {

    /**
     * Length of the compressed byte array.
     */
    public final int compressedLength;
    /**
     * Length of the uncompressed byte array.
     */
    public final int uncompressedLength;

    /**
     * Constructor.
     *
     * @param compressedLength Length of the compressed byte array
     * @param uncompressedLength Length of the uncompressed byte array
     */
    public IndexedGzCoordinates(
            int compressedLength,
            int uncompressedLength
    ) {

        this.compressedLength = compressedLength;
        this.uncompressedLength = uncompressedLength;

    }

}
