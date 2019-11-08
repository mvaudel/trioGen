package no.uib.triogen.io.genotypes;

/**
 * The variant iterator iterates through  all the variants of a file.
 *
 * @author Marc Vaudel
 */
public interface VariantIterator extends AutoCloseable {
    
    /**
     * Read the next variant. Returns null if there is no next variant.
     * 
     * @return a GenotypesProvider for the next variant
     */
    public GenotypesProvider next();

}
