package no.uib.triogen.processing.linear_model;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.getIndexFile;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.covariates.CovariatesHandler;
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
     * The allele inheritance map.
     */
    private final HashMap<Integer, char[]> inheritanceMap;
    /**
     * The default ploidy for mothers.
     */
    private final int defaultMotherPloidy;
    /**
     * The default ploidy for fathers.
     */
    private final int defaultFatherPloidy;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * The allele frequency threshold.
     */
    private final double alleleFrequencyThreshold;
    /**
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The names of the phenotypes to use.
     */
    private final String[] phenoNames;
    /**
     * Names of the covariates to use for all the phenotypes.
     */
    public final String[] covariatesGeneral;
    /**
     * Map of the covariates to use for specific phenotypes.
     */
    public final HashMap<String, TreeSet<String>> covariatesSpecific;
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
     * @param inheritanceMap The inheritance map for the given file.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     * @param variantList The variants to process.
     * @param mafThreshold The maf threshold.
     * @param childToParentMap The map of trios.
     * @param phenotypesFile The file containing the phenotypes.
     * @param phenoNames The names of the phenotypes to use.
     * @param covariatesGeneral The names of the general covariates to use for
     * all phenotypes.
     * @param covariatesSpecific The names of the covariates specific to
     * phenotypes.
     * @param models The models to use.
     * @param destinationFile The file to export the result to.
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public LinearModelComputer(
            File genotypesFile,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy,
            VariantList variantList,
            double mafThreshold,
            ChildToParentMap childToParentMap,
            File phenotypesFile,
            String[] phenoNames,
            String[] covariatesGeneral,
            HashMap<String, TreeSet<String>> covariatesSpecific,
            Model[] models,
            File destinationFile,
            int nVariants,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;
        this.variantList = variantList;
        this.alleleFrequencyThreshold = mafThreshold;
        this.childToParentMap = childToParentMap;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.covariatesGeneral = covariatesGeneral;
        this.covariatesSpecific = covariatesSpecific;
        this.models = models;
        this.destinationFile = destinationFile;
        this.nVariants = nVariants;
        this.logger = logger;

    }

    /**
     * Runs the linear association.
     *
     * @param timeOutDays The time out time in days.
     *
     * @throws InterruptedException Exception thrown if the process was.
     * interrupted
     * @throws TimeoutException Exception thrown if the process timed out.
     * @throws IOException Exception thrown if an i/o error occurred.
     */
    public void run(
            int timeOutDays
    ) throws InterruptedException, TimeoutException, IOException {

        logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        BgenIndex bgenIndex = BgenIndex.getBgenIndex(genotypesFile);
        BgenFileReader bgenFileReader  = new BgenFileReader(
                genotypesFile, 
                bgenIndex, 
                inheritanceMap, 
                defaultMotherPloidy, 
                defaultFatherPloidy
        );

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing " + genotypesFile + " done (" + duration + " seconds)");

        logger.logMessage("Importing " + phenoNames.length + " phenotyes from " + phenotypesFile.getAbsolutePath() + " for association with " + genotypesFile.getName());

        start = Instant.now().getEpochSecond();

        HashSet<String> allPhenoColumns = Stream.concat(
                Arrays.stream(phenoNames),
                Stream.concat(
                        Arrays.stream(covariatesGeneral),
                        Stream.concat(
                                covariatesSpecific.keySet().stream(),
                                covariatesSpecific.values().stream()
                                        .flatMap(
                                                subCovariates -> subCovariates.stream()
                                        )
                        )
                )
        ).collect(
                Collectors.toCollection(HashSet::new)
        );

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(
                phenotypesFile,
                childToParentMap.children,
                allPhenoColumns
        );

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (" + phenoNames.length + " phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)");

        logger.logMessage("Adjusting for covariates for association with " + genotypesFile.getName());

        start = Instant.now().getEpochSecond();

        CovariatesHandler covariatesHandler = new CovariatesHandler(
                phenoNames,
                phenotypesHandler,
                covariatesGeneral,
                covariatesSpecific
        );

        phenotypesHandler.adjustForCovariates(phenoNames, covariatesHandler);

        covariatesHandler.trim();

        end = Instant.now().getEpochSecond();
        duration = end - start;

        Arrays.stream(phenoNames)
                .sorted()
                .forEach(
                        phenoName -> logger.logMessage("    " + phenoName + ": " + covariatesHandler.originalIndexMap.get(phenoName).length + " phenotypes, " + covariatesHandler.covariatesMap.get(phenoName).length + " covariates (including intercept) representing " + covariatesHandler.rankMap.get(phenoName) + " dimensions.")
                );

        logger.logMessage("Done (Adjusted for covariates in " + duration + " seconds)");

        phenotypesHandler.sanityCheck();

        String nVariantsText = variantList == null ? "" : ", " + variantList.variantId.length + " variants";

        logger.logMessage("Linear association (geno: " + genotypesFile.getAbsolutePath() + nVariantsText + ", pheno: " + phenotypesFile.getAbsolutePath() + " " + phenoNames.length + " phenotypes)");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = new VariantIterator(
                bgenIndex,
                logger, 
                "Linear association in " + genotypesFile.getAbsolutePath(),
                variantList == null
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
                "rsId",
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
                        "rsId",
                        "testedAllele",
                        "otherAllele",
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
                                    bgenIndex,
                                    bgenFileReader,
                                    variantList,
                                    alleleFrequencyThreshold,
                                    childToParentMap,
                                    models,
                                    phenotypesHandler,
                                    covariatesHandler,
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

        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Done (Linear model for " + genotypesFile.getName() + ", " + bgenIndex.variantInformationArray.length + " variants and " + phenoNames.length + " phenoptyes processed in " + duration + " seconds)");

    }
}
