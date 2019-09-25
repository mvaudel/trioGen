package no.uib.triogen.io.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class iterates all the variants of a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfIterator {
    
    /**
     * The vcf file reader to use for iteration.
     */
    private final VCFFileReader vcfFileReader;
    
    private final CloseableIterator<VariantContext> iterator;
    /**
     * Mutex to query the reader.
     */
    private final SimpleSemaphore mutex = new SimpleSemaphore(1);

    /**
     * Constructor.
     * 
     * @param vcfFileReader the vcf file reader to use for iteration
     */
    public VcfIterator(VCFFileReader vcfFileReader) {
        
        this.vcfFileReader = vcfFileReader;
        iterator = vcfFileReader.iterator();
        
    }
    
    /**
     * Returns the next variant context.
     * 
     * @return the next variant context
     */
    public VariantContext next() {
        
        mutex.acquire();
        
        VariantContext variantContext = iterator.next();
        
        mutex.release();
        
        return variantContext;
        
    }
    
    /**
     * Closes the iterator
     */
    public void close() {
        
        iterator.close();
        
    }
    
}
