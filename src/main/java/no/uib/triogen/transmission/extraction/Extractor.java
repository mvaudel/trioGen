package no.uib.triogen.transmission.extraction;

import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
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
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;

    /**
     * Constructors.
     *
     * @param vcfFile the vcf file to use
     * @param childToParentMap the map of trios
     */
    public Extractor(
            File vcfFile,
            ChildToParentMap childToParentMap
    ) {

        this.vcfFile = vcfFile;
        this.childToParentMap = childToParentMap;

    }

    /**
     * Runs the extraction
     *
     * @param nThreads the number of threads to use
     * @param timeOutDays the time out time in days
     *
     * @throws InterruptedException exception thrown if the process was
     * interrupted
     * @throws TimeoutException exception thrown if the process timed out
     */
    public void run(int nThreads, int timeOutDays) throws InterruptedException, TimeoutException {

        VcfIterator iterator = new VcfIterator(vcfFile);

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nThreads);

            IntStream.range(0, nThreads)
                    .mapToObj(
                            i -> new ExtractorRunnable(iterator)
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
    }

}
