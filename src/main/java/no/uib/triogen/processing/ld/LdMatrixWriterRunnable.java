package no.uib.triogen.processing.ld;

import no.uib.triogen.io.genotypes.iterators.BufferedGenotypesIterator;
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
    private final BufferedGenotypesIterator iterator;
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
     * @param iterator The variant iterator.
     * @param logger The logger.
     */
    public LdMatrixWriterRunnable(
            BufferedGenotypesIterator iterator,
            Logger logger
    ) {
        
        this.iterator = iterator;
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
