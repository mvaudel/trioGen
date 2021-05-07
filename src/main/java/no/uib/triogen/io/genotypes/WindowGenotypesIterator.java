package no.uib.triogen.io.genotypes;

/**
 * Iterator for genotypes with a window around the variant being iterated.
 *
 * Marc Vaudel
 */
public interface WindowGenotypesIterator {

    /**
     * Returns a genotypes provider for the next variant in the iterator. Null if iteration is finished.
     * 
     * @return A genotypes provider for the next variant in the iterator. Null if iteration is finished.
     */
    public GenotypesProvider next();

    /**
     * Returns an iterator over the genotypes of the given contig in the given bp
     * range.
     *
     * @param contig The contig.
     * @param startBp The bp range start (inclusive).
     * @param endBp The bp range end (inclusive).
     *
     * @return An array of the genotypes.
     */
    public VariantIterator getGenotypesInRange(
            String contig,
            int startBp,
            int endBp
    );

    /**
     * Release a lower-bound for the section of the file to keep in memory.
     *
     * @param contig The contig.
     * @param minBp The base pair coordinate to release.
     */
    public void releaseMinBp(
            String contig,
            int minBp
    );
    
}