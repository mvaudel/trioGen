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
import no.uib.triogen.log.Logger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.Model;
import no.uib.triogen.model.geno.VariantList;
import no.uib.triogen.model.pheno.PhenotypesHandler;

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
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included.
     */
    private final double mafThreshold;
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
    private final Logger logger;

    /**
     * Constructor.
     *
     * @param genotypesFile the file containing the genotypes
     * @param genotypesFileType the type of genotypes file
     * @param variantList the variants to process
     * @param mafThreshold the maf threshold
     * @param childToParentMap the map of trios
     * @param phenotypesFile the file containing the phenotypes
     * @param phenoNames the names of the phenotypes to use
     * @param covariates the names of the covariates to use
     * @param models the models to use
     * @param destinationFile the file to export the result to
     * @param nVariants the number of variants to process in parallel
     * @param logger the logger
     */
    public LinearModelComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            VariantList variantList,
            double mafThreshold,
            ChildToParentMap childToParentMap,
            File phenotypesFile,
            String[] phenoNames,
            String[] covariates,
            Model[] models,
            File destinationFile,
            int nVariants,
            Logger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.variantList = variantList;
        this.mafThreshold = mafThreshold;
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
     * @param timeOutDays the time out time in days
     * @param test in test mode only a few variants will be processed
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
                variantList
        );
        IndexedGzWriter outputWriter = new IndexedGzWriter(
                destinationFile
        );

        File indexFile = getIndexFile(destinationFile);
        SimpleFileWriter index = new SimpleFileWriter(indexFile, true);
        index.writeLine(
                "VariantId",
                "Phenotype",
                "CompressedLength",
                "UncompressedLength"
        );

        IndexedGzCoordinates coordinates = outputWriter.append("# TrioGen version: " + TrioGen.getVersion() + IoUtils.lineSeparator);
        index.writeLine(
                "Header",
                "Header",
                Integer.toString(coordinates.compressedLength),
                Integer.toString(coordinates.uncompressedLength)
        );

        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(
                String.join(IoUtils.separator,
                        "phenotype",
                        "variantId",
                        "n",
                        "nAlt",
                        "nH"
                )
        );
        for (Model model : models) {

            model.getHeader(stringBuilder);

        }

        stringBuilder.append(IoUtils.lineSeparator);

        coordinates = outputWriter.append(stringBuilder.toString());
        index.writeLine(
                "Header",
                "Header",
                Integer.toString(coordinates.compressedLength),
                Integer.toString(coordinates.uncompressedLength)
        );

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LinearModelRunnable(
                                    iterator,
                                    mafThreshold,
                                    childToParentMap,
                                    models,
                                    phenotypesHandler,
                                    outputWriter,
                                    index,
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
