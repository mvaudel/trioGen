package no.uib.triogen.processing.ld;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.iterators.BufferedGenotypesIterator;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.Logger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.InputType;
import no.uib.triogen.model.geno.VariantIndex;
import no.uib.triogen.model.geno.VariantList;

/**
 * This class iterates through genotypes and writes a matrix of ld between
 * markers computed using sliding widows of the specified range.
 *
 * @author Marc Vaudel
 */
public class LdMatrixComputer {

    /**
     * Loading factor used to make sure that the buffer contains the ld range. A loading factor of 2 for a sliding window of 10 bp results in buffering 20 pb.
     */
    public static final double LOADING_FACTOR = 2.0;
    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The type of genotype file.
     */
    private final GenotypesFileType genotypesFileType;
    /**
     * The child ids of the trios to include.
     */
    private final String[] childIds;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The bp distance used to compute the ld in sliding windows. A max distance of 10 bp means a sliding window of 20 bp.
     */
    private final int maxDistance;
    /**
     * The genotype input type.
     */
    private final InputType inputType;
    /**
     * Index for the variants.
     */
    private final VariantIndex variantIndex = new VariantIndex();
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
     * @param childIds The child ids of the trios to include.
     * @param childToParentMap The map of trios.
     * @param variantList The variants to process.
     * @param destinationFile File to write to.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param inputType The genotype input type.
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public LdMatrixComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            String[] childIds,
            ChildToParentMap childToParentMap,
            VariantList variantList,
            File destinationFile,
            int maxDistance,
            InputType inputType,
            int nVariants,
            Logger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.childIds = childIds;
        this.childToParentMap = childToParentMap;
        this.variantList = variantList;
        this.destinationFile = destinationFile;
        this.maxDistance = maxDistance;
        this.inputType = inputType;
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
        
        BufferedGenotypesIterator bufferedIterator = new BufferedGenotypesIterator(
                iterator, 
                (int) LOADING_FACTOR * maxDistance, 
                (int) LOADING_FACTOR * maxDistance
        );
        
        LdMatrixWriter writer = new LdMatrixWriter(
                variantIndex,
                destinationFile
        );

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LdMatrixComputerRunnable(
                                    writer,
                                    bufferedIterator, 
                                    childIds, 
                                    childToParentMap, 
                                    maxDistance, 
                                    inputType,
                                    variantIndex, 
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
