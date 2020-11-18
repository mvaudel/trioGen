package no.uib.triogen.processing.ld;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.iterators.window.BufferedGenotypesIterator;
import no.uib.triogen.io.genotypes.iterators.window.TargetGenotypesIterator;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantIndex;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class iterates through genotypes and writes a matrix of ld between
 * markers computed using sliding widows of the specified range.
 *
 * @author Marc Vaudel
 */
public class LdMatrixComputer {

    /**
     * The loading factor can be used to reduce the frequency of buffering and
     * cache clean-up. With a loading factor of two, the buffer will be filled
     * twice what is needed.
     */
    private final double downstreamLoadingFactor;
    /**
     * The loading factor can be used to reduce the frequency of buffering and
     * cache clean-up. With a loading factor of two, the buffer will be filled
     * twice what is needed.
     */
    private final double upstreamLoadingFactor;
    /**
     * Loading factor used to make sure that the buffer contains the ld range. A
     * loading factor of 2 for a sliding window of 10 bp results in buffering 20
     * pb.
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
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The maximal bp distance used to compute the ld between variants. A max distance
     * of 10 bp means a window of 20 bp.
     */
    private final int maxDistance;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included (inclusive).
     */
    private final double mafThreshold;
    /**
     * Boolean indicating whether hard calls should be used.
     */
    private final boolean hardCalls;
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
     * @param genotypesFileType The type of genotypes file.
     * @param variantList The variants to process.
     * @param childToParentMap The map of trios.
     * @param destinationStem The stem of the path of the file where to write
     * the results.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param minR2 The minimal ld r2 to report (inclusive).
     * @param mafThreshold The maf threshold. maf is computed in parents and
     * values lower than threshold are not included (inclusive).
     * @param hardCalls Boolean indicating whether hard calls should be used
     * instead of dosages.
     * @param nVariants The number of variants to process in parallel.
     * @param downstreamLoadingFactor The loading factor to use when trimming
     * the buffer.
     * @param upstreamLoadingFactor The loading factor to use when buffering.
     * @param logger The logger.
     */
    public LdMatrixComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            VariantList variantList,
            ChildToParentMap childToParentMap,
            String destinationStem,
            int maxDistance,
            double minR2,
            double mafThreshold,
            boolean hardCalls,
            int nVariants,
            double downstreamLoadingFactor,
            double upstreamLoadingFactor,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.destinationStem = destinationStem;
        this.maxDistance = maxDistance;
        this.minR2 = minR2;
        this.mafThreshold = mafThreshold;
        this.hardCalls = hardCalls;
        this.nVariants = nVariants;
        this.downstreamLoadingFactor = downstreamLoadingFactor;
        this.upstreamLoadingFactor = upstreamLoadingFactor;
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

        long start = Instant.now().getEpochSecond();
        
        String readout = hardCalls ? "genotypes" : "dosages"; 

        if (variantList == null) {

            logger.logMessage(genotypesFile.getName() + " LD extraction on " + readout + " using " + nVariants + " threads.");

            VariantIterator iterator = GenotypesFileType.getVariantIterator(
                    genotypesFile,
                    genotypesFileType,
                    logger
            );

            BufferedGenotypesIterator bufferedIterator = new BufferedGenotypesIterator(
                    iterator,
                    childToParentMap,
                    (int) LOADING_FACTOR * maxDistance,
                    (int) LOADING_FACTOR * maxDistance,
                    mafThreshold,
                    nVariants,
                    downstreamLoadingFactor,
                    upstreamLoadingFactor,
                    true
            );

            File destinationFile = new File(destinationStem + ".tld");

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
                                        childToParentMap,
                                        maxDistance,
                                        minR2,
                                        hardCalls,
                                        variantIndex,
                                        testIteration,
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

                writer.close();
                iterator.close();

            }

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            logger.logMessage(genotypesFile.getName() + " Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)");

        } else {

            int nVariants = variantList.variantId.length;

            logger.logMessage(genotypesFile.getName() + " LD extraction (" + nVariants + " target variants) on " + readout + " using " + nVariants + " threads.");

            ArrayList<LdMatrixWriter> writers = new ArrayList<>();

            try {

                ExecutorService pool = Executors.newFixedThreadPool(nVariants);

                for (int i = 0; i < nVariants; i++) {

                    String variantId = variantList.variantId[i];

                    TargetGenotypesIterator iterator = new TargetGenotypesIterator(
                            genotypesFile,
                            variantId,
                            variantList.chromosome[i],
                            variantList.start[i],
                            variantList.end[i],
                            maxDistance,
                            maxDistance,
                            mafThreshold,
                            childToParentMap
                    );

                    File destinationFile = new File(destinationStem + "_" + variantId + ".tld");

                    LdMatrixWriter writer = new LdMatrixWriter(
                            variantIndex,
                            destinationFile
                    );
                    writers.add(writer);

                    LdMatrixComputerRunnable runnable = new LdMatrixComputerRunnable(
                            writer,
                            iterator,
                            childToParentMap,
                            maxDistance,
                            minR2,
                            hardCalls,
                            variantIndex,
                            testIteration,
                            logger
                    );

                    pool.submit(runnable);

                }

                pool.shutdown();

                if (!pool.awaitTermination(timeOutDays, TimeUnit.DAYS)) {

                    throw new TimeoutException("Analysis timed out (time out: " + timeOutDays + " days)");

                }

            } finally {

                writers.stream()
                        .forEach(
                                iterator -> {
                                    try {
                                        iterator.close();
                                    } catch (Throwable t) {
                                        t.printStackTrace();
                                    }
                                }
                        );

            }

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            logger.logMessage(genotypesFile.getName() + " Done (" + nVariants + " variants processed in " + duration + " seconds)");

        }
    }
}
