package no.uib.triogen.io.genotypes;

/**
 * The variant iterator iterates through  all the variants of a file.
 *
 * @author Marc Vaudel
 */
public interface VariantIterator {
    
    /**
     * Read the next variant. Returns null if there is no next variant.
     * 
     * @return a GenotypesProvider for the next variant
     */
    public GenotypesProvider next();

    /**
     * Returns the number of variants read from the file.
     *
     * @return the number of variants read from the file
     */
    public int getnVariants();

    /**
     * Closes the iterator.
     */
    public void close();
    
    /**
     * Returns a boolean indicating whether the iteration is finished. 
     * 
     * @return A boolean indicating whether the iteration is finished.
     */
    public boolean isFinished();

}
