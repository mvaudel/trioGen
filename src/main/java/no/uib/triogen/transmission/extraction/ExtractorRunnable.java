package no.uib.triogen.transmission.extraction;

import htsjdk.variant.variantcontext.VariantContext;
import no.uib.triogen.io.vcf.VcfIterator;

/**
 * Runnable for the extraction of transmitted alleles.
 *
 * @author Marc Vaudel
 */
public class ExtractorRunnable implements Runnable {

    private final VcfIterator iterator;

    public ExtractorRunnable(
            VcfIterator iterator
    ) {

        this.iterator = iterator;

    }

    @Override
    public void run() {
        
        
    }
}
