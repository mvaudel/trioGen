package no.uib.triogen.transmission.extraction;

import java.io.File;
import java.time.Instant;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.vcf.VcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * This class iterates a vcf file and extracts transmitted alleles.
 *
 * @author Marc Vaudel
 */
public class Extractor {

    /**
     * The vcf file to use.
     */
    private final File vcfFile;
    /**
     * The file to export the result to.
     */
    private final File destinationFile;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The number of variants to process at the same time.
     */
    private final int nThreads = 4;

    /**
     * Constructors.
     *
     * @param vcfFile the vcf file to use
     * @param childToParentMap the map of trios
     * @param destinationFile the file to export the result to
     */
    public Extractor(
            File vcfFile,
            ChildToParentMap childToParentMap,
            File destinationFile
    ) {

        this.vcfFile = vcfFile;
        this.childToParentMap = childToParentMap;
        this.destinationFile = destinationFile;

    }

    /**
     * Runs the extraction
     *
     * @param timeOutDays the time out time in days
     * @param test in test mode only a few variants will be processed
     *
     * @throws InterruptedException exception thrown if the process was
     * interrupted
     * @throws TimeoutException exception thrown if the process timed out
     */
    public void run(
            int timeOutDays,
            boolean test
    ) throws InterruptedException, TimeoutException {

        System.out.println(
                Instant.now() + " - Starting processing of " + vcfFile.getAbsolutePath()
        );

        long start = Instant.now().getEpochSecond();

        VcfIterator iterator = new VcfIterator(
                vcfFile
        );
        SimpleFileWriter writer = new SimpleFileWriter(
                destinationFile,
                test
        );

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nThreads);

            IntStream.range(0, nThreads)
                    .mapToObj(
                            i -> new ExtractorRunnable(
                                    iterator,
                                    childToParentMap,
                                    writer
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

        System.out.println(
                Instant.now() + " - Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)"
        );
    }
}
