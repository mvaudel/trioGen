package no.uib.triogen.processing.ld;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.genotypes.bgen.VariantIterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantIndex;

/**
 * This class iterates through genotypes and writes a matrix of ld between
 * markers computed using sliding widows of the specified range.
 *
 * @author Marc Vaudel
 */
public class LdMatrixComputer {

    /**
     * The file containing the genotypes.
     */
    private final File genotypesFile;
    /**
     * The inheritance map.
     */
    private final HashMap<Integer, char[]> inheritanceMap;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The maximal bp distance used to compute the ld between variants. A max
     * distance of 10 bp means a window of 20 bp.
     */
    private final int maxDistance;
    /**
     * Index for the variants.
     */
    private final VariantIndex variantIndex = new VariantIndex();
    /**
     * The stem of the path of the file where to write the results.
     */
    private final String destinationStem;
    /**
     * The number of variants to process in parallel.
     */
    private final int nVariants;
    /**
     * The minimal ld r2 to report (inclusive).
     */
    private final double minR2;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * constructor.
     *
     * @param genotypesFile The file containing the genotypes.
     * @param inheritanceMap The inheritance map for the given file.
     * @param childToParentMap The map of trios.
     * @param destinationStem The stem of the path of the file where to write
     * the results.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param minR2 The minimal ld r2 to report (inclusive).
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public LdMatrixComputer(
            File genotypesFile,
            HashMap<Integer, char[]> inheritanceMap,
            ChildToParentMap childToParentMap,
            String destinationStem,
            int maxDistance,
            double minR2,
            int nVariants,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.inheritanceMap = inheritanceMap;
        this.childToParentMap = childToParentMap;
        this.destinationStem = destinationStem;
        this.maxDistance = maxDistance;
        this.minR2 = minR2;
        this.nVariants = nVariants;
        this.logger = logger;

    }

    /**
     * Runs the matrix export.
     *
     * @param timeOutDays The time out time in days.
     * @param testIteration In testIteration mode LD calculations will not be
     * computed.
     * @param test In test mode only a few variants will be processed.
     *
     * @throws InterruptedException exception thrown if the process was
     * interrupted
     * @throws TimeoutException exception thrown if the process timed out
     * @throws IOException exception thrown if an i/o error occurred
     */
    public void run(
            int timeOutDays,
            boolean testIteration,
            boolean test
    ) throws InterruptedException, TimeoutException, IOException {

        if (test) {

            logger.logMessage("*** TEST MODE ***");

        }

        logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        BgenIndex bgenIndex = BgenIndex.getBgenIndex(genotypesFile);
        BgenFileReader bgenFileReader = new BgenFileReader(genotypesFile, bgenIndex, inheritanceMap);

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing " + genotypesFile + " done (" + duration + " seconds)");

        logger.logMessage("LD extraction in " + genotypesFile.getName() + " using " + nVariants + " threads.");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = new VariantIterator(bgenIndex);

        P0Cache p0Cache = new P0Cache(nVariants);

        File destinationFile = new File(destinationStem + ".tld");

        try (
                LdMatrixWriter writer = new LdMatrixWriter(
                        variantIndex,
                        destinationFile
                )) {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LdMatrixComputerRunnable(
                                    writer,
                                    iterator,
                                    bgenIndex,
                                    bgenFileReader,
                                    childToParentMap,
                                    maxDistance,
                                    minR2,
                                    variantIndex,
                                    p0Cache,
                                    i,
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

        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage(genotypesFile.getName() + " Done (" + bgenIndex.variantInformationArray.length + " variants processed in " + duration + " seconds)");

    }
}
