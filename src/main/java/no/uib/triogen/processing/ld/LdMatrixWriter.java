package no.uib.triogen.processing.ld;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
import no.uib.triogen.log.Logger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.VariantList;
import no.uib.triogen.processing.association.linear_model.LinearModelRunnable;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class iterates through genotypes and writes a matrix of ld between
 * markers.
 *
 * @author Marc Vaudel
 */
public class LdMatrixWriter {

    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The type of genotype file.
     */
    private final GenotypesFileType genotypesFileType;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The maximal number of bp to allow between variants.
     */
    private final int maxDistance;
    /**
     * The file where to write the results.
     */
    private final File destinationFile;
    /**
     * The number of variants to process in parallel.
     */
    private final int nVariants;
    /**
     * The logger.
     */
    private final Logger logger;

    /**
     * constructor.
     *
     * @param genotypesFile The file containing the genotypes.
     * @param genotypesFileType The type of genotypes file.
     * @param childToParentMap The map of trios.
     * @param variantList The variants to process.
     * @param destinationFile File to write to.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public LdMatrixWriter(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            ChildToParentMap childToParentMap,
            VariantList variantList,
            File destinationFile,
            int maxDistance,
            int nVariants,
            Logger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.childToParentMap = childToParentMap;
        this.variantList = variantList;
        this.destinationFile = destinationFile;
        this.maxDistance = maxDistance;
        this.nVariants = nVariants;
        this.logger = logger;

    }

    /**
     * Runs the matrix export.
     *
     * @param timeOutDays The time out time in days.
     * @param test In test mode only a few variants will be processed.
     *
     * @throws InterruptedException exception thrown if the process was
     * interrupted
     * @throws TimeoutException exception thrown if the process timed out
     * @throws IOException exception thrown if an i/o error occurred
     */
    public void run(
            int timeOutDays,
            boolean test
    ) throws InterruptedException, TimeoutException, IOException {

        if (test) {

            logger.logMessage("*** TEST MODE ***");

        }

        String nVariantsText = variantList == null ? "" : ", " + variantList.variantId.length + " variants";

        logger.logMessage("LD extraction (geno: " + genotypesFile.getAbsolutePath() + nVariantsText + ")");

        long start = Instant.now().getEpochSecond();

        VariantIterator iterator = GenotypesFileType.getVariantIterator(
                genotypesFile,
                genotypesFileType,
                variantList
        );
        
        

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LdMatrixWriterRunnable(
                                    iterator,
                                    mutex,
                                    currentBp,
                                    variantIndex,
                                    buffer,
                                    logger
                            )
                    )
                    .forEach(
                            worker -> pool.submit(worker)
                    );

            pool.shutdown();

            if (!pool.awaitTermination(timeOutDays, TimeUnit.DAYS)) {

                throw new TimeoutException("Analysis timed out (time out: " + timeOutDays + " days)");

            }

        } finally {

            iterator.close();

        }

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)");

    }

}
