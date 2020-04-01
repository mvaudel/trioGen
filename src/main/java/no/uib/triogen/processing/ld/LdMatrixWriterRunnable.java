package no.uib.triogen.processing.ld;

import no.uib.triogen.io.genotypes.iterators.BufferredGenotypesIterator;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
import no.uib.triogen.log.Logger;

/**
 * Runnable for the LD matrix writer.
 *
 * @author Marc Vaudel
 */
public class LdMatrixWriterRunnable implements Runnable {

    /**
     * The buffer.
     */
    private final BufferredGenotypesIterator buffer;
    /**
     * The logger.
     */
    private final Logger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;
    
    /**
     * Constructor.
     * 
     * @param buffer The variant buffer.
     * @param logger The logger.
     */
    public LdMatrixWriterRunnable(
            BufferredGenotypesIterator buffer,
            Logger logger
    ) {
        
        this.buffer = buffer;
        this.logger = logger;
        
    }

    @Override
    public void run() {

        try {

            

        } catch (Throwable t) {

            canceled = true;

            logger.logError(
                    Arrays.stream(t.getStackTrace())
                            .map(
                                    element -> element.toString()
                            )
                            .collect(Collectors.joining(" "))
            );

            t.printStackTrace();

        }
    }
}
