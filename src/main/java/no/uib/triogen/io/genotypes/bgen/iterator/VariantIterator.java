package no.uib.triogen.io.genotypes.bgen.iterator;

import java.time.Instant;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Iterator for the variants of a bgen file.
 *
 * @author Marc Vaudel
 */
public class VariantIterator {

    /**
     * The number of variants to iterate before showing progress.
     */
    private static final int nProgress = 100000;
    /**
     * The index of the bgen file to iterate.
     */
    private final BgenIndex bgenIndex;
    /**
     * The position to start the iteration. Ignored if -1.
     */
    private final int start;
    /**
     * The position to end the iteration. Ignored if -1.
     */
    private final int end;
    /**
     * A semaphore to synchronize threads.
     */
    private final SimpleSemaphore simpleSemaphore = new SimpleSemaphore(1);
    /**
     * The current index of the iteration.
     */
    private int currentVariantIndex = 0;
    /**
     * The logger to use. Ignored if null.
     */
    private final SimpleCliLogger logger;
    /**
     * The prefix to use for the log.
     */
    private final String logPrefix;
    /**
     * Instant where the iteration began.
     */
    private long startInstant = -1l;
    /**
     * Boolean indicating whether estimates for time of completion should be
     * given.
     */
    private final boolean estimateTime;

    /**
     * Constructor.
     *
     * @param bgenIndex The index of the bgen file to iterate.
     * @param start The position to start the iteration. Ignored if -1.
     * @param end The position to end the iteration. Ignored if -1.
     * @param logger The logger to use. Ignored if null.
     * @param logPrefix The prefix to use for the log.
     * @param estimateTime Boolean indicating whether estimates for time of
     * completion should be given.
     */
    public VariantIterator(
            BgenIndex bgenIndex,
            int start,
            int end,
            SimpleCliLogger logger,
            String logPrefix,
            boolean estimateTime
    ) {

        this.bgenIndex = bgenIndex;
        this.start = start;
        this.end = end;
        this.logger = logger;
        this.logPrefix = logPrefix;
        this.estimateTime = estimateTime;

    }

    /**
     * Constructor.
     *
     * @param bgenIndex The index of the bgen file to iterate.
     * @param start The position to start the iteration. Ignored if -1.
     * @param end The position to end the iteration. Ignored if -1.
     */
    public VariantIterator(
            BgenIndex bgenIndex,
            int start,
            int end
    ) {

        this(bgenIndex, start, end, null, null, false);

    }

    /**
     * Constructor.
     *
     * @param bgenIndex The index of the bgen file to iterate.
     * @param logger The logger to use. Ignored if null.
     * @param logPrefix The prefix to use for the log.
     * @param estimateTime Boolean indicating whether estimates for time of
     * completion should be given.
     */
    public VariantIterator(
            BgenIndex bgenIndex,
            SimpleCliLogger logger,
            String logPrefix,
            boolean estimateTime
    ) {

        this(bgenIndex, -1, -1, logger, logPrefix, estimateTime);

    }

    /**
     * Constructor.
     *
     * @param bgenIndex The index of the bgen file to iterate.
     */
    public VariantIterator(
            BgenIndex bgenIndex
    ) {

        this(bgenIndex, -1, -1, null, null, false);

    }

    /**
     * Returns the next position.
     *
     * @return The next position.
     */
    public Integer next() {

        simpleSemaphore.acquire();

        if (logger != null && currentVariantIndex % nProgress == 0) {

            double progress = ((double) (Math.round(1000.0 * currentVariantIndex) / bgenIndex.variantInformationArray.length)) / 10;

            if (estimateTime) {

                if (startInstant == -1l) {

                    startInstant = Instant.now().getEpochSecond();

                } else {

                    double elapsedTimeHours = ((double) Instant.now().getEpochSecond() - startInstant) / 3600;
                    double elapsedTimeHoursRounded = Math.round(elapsedTimeHours * 10) / 10;

                    double timeRemaining = ((double) bgenIndex.variantInformationArray.length) / currentVariantIndex * elapsedTimeHours;
                    double timeRemainingRounded = Math.round(timeRemaining * 10) / 10;

                    logger.logMessage(logPrefix + "    " + currentVariantIndex + " processed of " + bgenIndex.variantInformationArray.length + " (" + progress + "% in " + elapsedTimeHoursRounded + " hours, approx " + timeRemainingRounded + " hours remaining)");

                }

            } else if (currentVariantIndex > 0) {

                logger.logMessage(logPrefix + "    " + currentVariantIndex + " processed of " + bgenIndex.variantInformationArray.length + " (" + progress + "%)");

            }
        }

        if ((start != -1 || end != -1) && currentVariantIndex < bgenIndex.variantInformationArray.length) {

            VariantInformation variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

            while (start != -1 && variantInformation.position < start || end != -1 && variantInformation.position > end) {

                currentVariantIndex++;

                if (currentVariantIndex != bgenIndex.variantInformationArray.length) {

                    variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

                } else {

                    break;

                }
            }
        }

        int index = currentVariantIndex;

        currentVariantIndex++;

        simpleSemaphore.release();

        return index < bgenIndex.variantInformationArray.length ? index : null;

    }

}
