package no.uib.triogen.processing.association.linear_model;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.getIndexFile;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class runs a linear model on the given phenotypes.
 *
 * @author Marc Vaudel
 */
public class LinearModelComputer {

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
     * Variants with a distance to a target variant lower than maxDistance will be included in the computation.
     */
    private final int maxDistance;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included.
     */
    private final double mafThreshold;
    /**
     * If true, dosages will be used where possible, hard calls otherwise.
     */
    private final boolean useDosages;
    /**
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The names of the phenotypes to use.
     */
    private final String[] phenoNames;
    /**
     * The covariates to use.
     */
    private final String[] covariates;
    /**
     * The models to use.
     */
    private final Model[] models;
    /**
     * The file to export the result to.
     */
    private final File destinationFile;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The number of variants to process in parallel.
     */
    private final int nVariants;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * Constructor.
     *
     * @param genotypesFile The file containing the genotypes.
     * @param genotypesFileType The type of genotypes file.
     * @param variantList The variants to process.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param mafThreshold The maf threshold.
     * @param useDosages If true, dosages will be used when possible, hard calls
     * otherwise.
     * @param childToParentMap The map of trios.
     * @param phenotypesFile The file containing the phenotypes.
     * @param phenoNames The names of the phenotypes to use.
     * @param covariates The names of the covariates to use.
     * @param models The models to use.
     * @param destinationFile The file to export the result to.
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public LinearModelComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            VariantList variantList,
            int maxDistance,
            double mafThreshold,
            boolean useDosages,
            ChildToParentMap childToParentMap,
            File phenotypesFile,
            String[] phenoNames,
            String[] covariates,
            Model[] models,
            File destinationFile,
            int nVariants,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.variantList = variantList;
        this.maxDistance = maxDistance;
        this.mafThreshold = mafThreshold;
        this.useDosages = useDosages;
        this.childToParentMap = childToParentMap;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.covariates = covariates;
        this.models = models;
        this.destinationFile = destinationFile;
        this.nVariants = nVariants;
        this.logger = logger;

    }

    /**
     * Runs the linear association.
     *
     * @param timeOutDays The time out time in days.
     * @param test In test mode only a few variants will be processed.
     *
     * @throws InterruptedException Exception thrown if the process was.
     * interrupted
     * @throws TimeoutException Exception thrown if the process timed out.
     * @throws IOException Exception thrown if an i/o error occurred.
     */
    public void run(
            int timeOutDays,
            boolean test
    ) throws InterruptedException, TimeoutException, IOException {

        if (test) {

            logger.logMessage("*** TEST MODE ***");

        }

        if (covariates.length > 0) {

            logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath() + " and adjusting for " + covariates.length + " covariates");

        } else {

            logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath());

        }

        long start = Instant.now().getEpochSecond();

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                phenoNames,
                covariates
        );

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Done (" + phenoNames.length + " phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)");

        String nVariantsText = variantList == null ? "" : ", " + variantList.variantId.length + " variants";

        logger.logMessage("Linear association (geno: " + genotypesFile.getAbsolutePath() + nVariantsText + ", pheno: " + phenotypesFile.getAbsolutePath() + ")");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = GenotypesFileType.getVariantIterator(
                genotypesFile,
                genotypesFileType,
                variantList,
                maxDistance
        );
        IndexedGzWriter outputWriter = new IndexedGzWriter(
                destinationFile
        );

        File indexFile = getIndexFile(destinationFile);
        SimpleFileWriter index = new SimpleFileWriter(indexFile, true);
        index.writeLine(
                "contig",
                "position",
                "variantId",
                "phenotype",
                "compressedLength",
                "uncompressedLength"
        );

        IndexedGzCoordinates coordinates = outputWriter.append("# TrioGen version: " + TrioGen.getVersion() + IoUtils.LINE_SEPARATOR);
        index.writeLine(
                "Header",
                "Header",
                "Header",
                "Comment",
                Integer.toString(coordinates.compressedLength),
                Integer.toString(coordinates.uncompressedLength)
        );

        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(
                String.join(IoUtils.SEPARATOR,
                        "phenotype",
                        "contig",
                        "position",
                        "variantId",
                        "refAllele",
                        "altAllele",
                        "typed",
                        "n",
                        "nAlt",
                        "nH",
                        "mendelianError"
                )
        );
        for (Model model : models) {

            model.getHeader(stringBuilder);

        }

        stringBuilder.append(IoUtils.LINE_SEPARATOR);

        coordinates = outputWriter.append(stringBuilder.toString());
        index.writeLine(
                "Header",
                "Header",
                "Header",
                "Header",
                Integer.toString(coordinates.compressedLength),
                Integer.toString(coordinates.uncompressedLength)
        );

        SimpleSemaphore gzIndexMutex = new SimpleSemaphore(1);

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LinearModelRunnable(
                                    iterator,
                                    variantList,
                                    mafThreshold,
                                    useDosages,
                                    childToParentMap,
                                    models,
                                    phenotypesHandler,
                                    outputWriter,
                                    index,
                                    gzIndexMutex,
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

            outputWriter.close();
            index.close();
            iterator.close();

        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)");

    }
}
