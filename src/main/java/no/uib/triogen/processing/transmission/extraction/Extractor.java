package no.uib.triogen.processing.transmission.extraction;

import java.io.File;
import java.time.Instant;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class iterates a vcf file and extracts transmitted alleles.
 *
 * @author Marc Vaudel
 */
public class Extractor {

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
     * The stem of the files to export the result to.
     */
    private final String destinationStem;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The number of variants to process in parallel.
     */
    private final int nVariants;

    /**
     * Constructor.
     *
     * @param genotypesFile the file containing the genotypes
     * @param genotypesFileType the type of genotypes file
     * @param variantList the variants to process
     * @param childToParentMap the map of trios
     * @param destinationStem the stem of the files to export the result to
     * @param nVariants the number of variants to process in parallel
     */
    public Extractor(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            VariantList variantList,
            ChildToParentMap childToParentMap,
            String destinationStem,
            int nVariants
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.destinationStem = destinationStem;
        this.nVariants = nVariants;

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

        if (test) {

            System.out.println("*** TEST MODE ***");

        }

        System.out.println(Instant.now() + " - Extraction of H from " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        VariantIterator iterator = GenotypesFileType.getVariantIterator(
                genotypesFile,
                genotypesFileType,
                variantList,
                false
        );
        SimpleFileWriter h1Writer = new SimpleFileWriter(
                new File(destinationStem + "_h1.gz"),
                true
        );
        SimpleFileWriter h2Writer = new SimpleFileWriter(
                new File(destinationStem + "_h2.gz"),
                true
        );
        SimpleFileWriter h3Writer = new SimpleFileWriter(
                new File(destinationStem + "_h3.gz"),
                true
        );
        SimpleFileWriter h4Writer = new SimpleFileWriter(
                new File(destinationStem + "_h4.gz"),
                true
        );

        h1Writer.writeLine("# TrioGen version: " + TrioGen.getVersion());
        h2Writer.writeLine("# TrioGen version: " + TrioGen.getVersion());
        h3Writer.writeLine("# TrioGen version: " + TrioGen.getVersion());
        h4Writer.writeLine("# TrioGen version: " + TrioGen.getVersion());

        String header = String.join(IoUtils.SEPARATOR,
                "contig",
                "position",
                "variantId",
                "refAllele",
                "altAllele",
                "genotyped",
                String.join(
                        IoUtils.SEPARATOR,
                        childToParentMap.children
                )
        );
        h1Writer.writeLine(header);
        h2Writer.writeLine(header);
        h3Writer.writeLine(header);
        h4Writer.writeLine(header);

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new ExtractorRunnable(
                                    iterator,
                                    childToParentMap,
                                    h1Writer,
                                    h2Writer,
                                    h3Writer,
                                    h4Writer
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

            h1Writer.close();
            h2Writer.close();
            h3Writer.close();
            h4Writer.close();
            iterator.close();

        }

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        System.out.println(
                Instant.now() + " - Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)"
        );
    }
}
