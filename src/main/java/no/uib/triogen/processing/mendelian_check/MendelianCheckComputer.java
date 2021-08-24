package no.uib.triogen.processing.mendelian_check;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class iterates through genotypes and writes a matrix of ld between
 * markers computed using sliding widows of the specified range.
 *
 * @author Marc Vaudel
 */
public class MendelianCheckComputer {

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
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included (inclusive).
     */
    private final double alleleFrequencyThreshold;
    /**
     * The path of the file where to write the results.
     */
    private final File destinationFile;
    /**
     * The number of variants to process in parallel.
     */
    private final int nVariants;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * constructor.
     *
     * @param genotypesFile The file containing the genotypes.
     * @param inheritanceMap The inheritance map for the given file.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     * @param variantList The variants to process.
     * @param childToParentMap The map of trios.
     * @param destinationFile The path of the file where to write the results.
     * @param alleleFrequencyThreshold The allele frequency threshold.
     * @param nVariants The number of variants to process in parallel.
     * @param logger The logger.
     */
    public MendelianCheckComputer(
            File genotypesFile,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy,
            VariantList variantList,
            ChildToParentMap childToParentMap,
            File destinationFile,
            double alleleFrequencyThreshold,
            int nVariants,
            SimpleCliLogger logger
    ) {

        this.genotypesFile = genotypesFile;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;
        this.variantList = variantList;
        this.childToParentMap = childToParentMap;
        this.destinationFile = destinationFile;
        this.alleleFrequencyThreshold = alleleFrequencyThreshold;
        this.nVariants = nVariants;
        this.logger = logger;

    }

    /**
     * Runs the matrix export.
     *
     * @param timeOutDays The time out time in days.
     *
     * @throws InterruptedException exception thrown if the process was
     * interrupted
     * @throws TimeoutException exception thrown if the process timed out
     * @throws IOException exception thrown if an i/o error occurred
     */
    public void run(
            int timeOutDays
    ) throws InterruptedException, TimeoutException, IOException {

        logger.logMessage("Parsing " + genotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        BgenIndex bgenIndex = BgenIndex.getBgenIndex(genotypesFile);
        BgenFileReader bgenFileReader = new BgenFileReader(
                genotypesFile,
                bgenIndex,
                inheritanceMap,
                defaultMotherPloidy,
                defaultFatherPloidy
        );

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing " + genotypesFile + " done (" + duration + " seconds)");

        start = Instant.now().getEpochSecond();

        if (variantList == null) {

            logger.logMessage("Mendelian error check in " + genotypesFile.getName() + " using " + nVariants + " threads.");

            VariantIterator iterator = new VariantIterator(
                    bgenIndex,
                    logger,
                    "Mendelian error check in " + genotypesFile.getAbsolutePath(),
                    variantList == null
            );

            try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

                writeHeader(writer);

                ExecutorService pool = Executors.newFixedThreadPool(nVariants);

                IntStream.range(0, nVariants)
                        .mapToObj(
                                i -> new MendelianCheckRunnable(
                                        writer,
                                        iterator,
                                        bgenIndex,
                                        bgenFileReader,
                                        childToParentMap,
                                        alleleFrequencyThreshold,
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

    /**
     * Writes the header to the result file.
     *
     * @param writer Writer to the result file
     */
    private static void writeHeader(
            SimpleFileWriter writer
    ) {

        writer.writeLine("# Medelian error check using TrioGen v." + TrioGen.getVersion());

        writer.writeLine(
                "contig",
                "position",
                "variantId",
                "refAllele",
                "altAllele",
                "typed",
                "maf",
                "h2_1_obs",
                "h2_2_obs",
                "h4_1_obs",
                "h4_2_obs",
                "h2_1_exp",
                "h2_2_exp",
                "h4_1_exp",
                "h4_2_exp",
                "h2_1_p",
                "h2_2_p",
                "h4_1_p",
                "h4_2_p",
                "prevalence_before_check",
                "prevalence_after_check"
        );

    }
}
