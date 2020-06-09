package no.uib.triogen.io.genotypes.iterators;

import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Genotype iterator based on an array.
 *
 * @author Marc Vaudel
 */
public class SimpleGenotypeIterator implements VariantIterator {
    
    /**
     * The array to iterate.
     */
    private final GenotypesProvider[] array;
    /**
     * The current iteration index.
     */
    private int currentIndex = 0;
    /**
     * Semaphore for currentIndex.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);
    /**
     * Constructor.
     * 
     * @param array The array to iterate.
     */
    public SimpleGenotypeIterator(
            GenotypesProvider[] array
    ) {
        
        this.array = array;
        
    }

    @Override
    public GenotypesProvider next() {
        
        semaphore.acquire();
        
        if (currentIndex >= array.length) {
            
            return null;
            
        }
        
        GenotypesProvider result = array[currentIndex];
        
        currentIndex++;
        
        semaphore.release();
        
        return result;
        
    }

    @Override
    public int getnVariants() {
        
        return currentIndex;
        
    }

    @Override
    public void close() {
        
        // Nothing to do here.
        
    }

    @Override
    public boolean isFinished() {
        
        return currentIndex >= array.length;
        
    }
            

}
