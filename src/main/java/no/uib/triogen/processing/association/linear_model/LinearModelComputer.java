package no.uib.triogen.processing.association.linear_model;

import java.io.File;
import java.time.Instant;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.IntStream;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
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
     * The file containing the phenotypes.
     */
    private final File phenotypesFile;
    /**
     * The names of the phenotypes to use.
     */
    private final String[] phenoNames;
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
     * Constructor.
     *
     * @param genotypesFile the file containing the genotypes
     * @param genotypesFileType the type of genotypes file
     * @param childToParentMap the map of trios
     * @param phenotypesFile the file containing the phenotypes
     * @param phenoNames the names of the phenotypes to use
     * @param destinationFile the file to export the result to
     * @param nVariants the number of variants to process in parallel
     */
    public LinearModelComputer(
            File genotypesFile,
            GenotypesFileType genotypesFileType,
            ChildToParentMap childToParentMap,
            File phenotypesFile,
            String[] phenoNames,
            File destinationFile,
            int nVariants
    ) {

        this.genotypesFile = genotypesFile;
        this.genotypesFileType = genotypesFileType;
        this.childToParentMap = childToParentMap;
        this.phenotypesFile = phenotypesFile;
        this.phenoNames = phenoNames;
        this.destinationFile = destinationFile;
        this.nVariants = nVariants;

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
     */
    public void run(
            int timeOutDays,
            boolean test
    ) throws InterruptedException, TimeoutException {

        if (test) {

            System.out.println("*** TEST MODE ***");

        }

        System.out.println(Instant.now() + " - Parsing phenotyes from " + phenotypesFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        PhenotypesHandler phenotypesHandler = new PhenotypesHandler(phenotypesFile, childToParentMap.children, phenoNames);

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        System.out.println(
                Instant.now() + " - Done (" + phenoNames.length + " phenotypes for " + phenotypesHandler.nChildren + " children imported in " + duration + " seconds)"
        );

        System.out.println(Instant.now() + " - Linear association (geno: " + genotypesFile.getAbsolutePath() + ", pheno: " + phenotypesFile.getAbsolutePath() + ")");

        start = Instant.now().getEpochSecond();

        VariantIterator iterator = GenotypesFileType.getVariantIterator(genotypesFile, genotypesFileType);
        SimpleFileWriter outputWriter = new SimpleFileWriter(
                destinationFile,
                true
        );

            String line = String.join(
                    Utils.separator,
                    "phenotype",
                    "variantId",
                    "nH",
                    "n",
                    "h1_beta",
                    "h1_beta_se",
                    "h1_beta_p",
                    "h2_beta",
                    "h2_beta_se",
                    "h2_beta_p",
                    "h3_beta",
                    "h3_beta_se",
                    "h3_beta_p",
                    "h4_beta",
                    "h4_beta_se",
                    "h4_beta_p",
                    "cmf_h_p",
                    "cmf_child_beta",
                    "cmf_child_beta_se",
                    "cmf_child_beta_p",
                    "cmf_mother_beta",
                    "cmf_mother_beta_se",
                    "cmf_mother_beta_p",
                    "cmf_father_beta",
                    "cmf_father_beta_se",
                    "cmf_father_beta_p",
                    "cm_cmf_p",
                    "cm_child_beta",
                    "cm_child_beta_se",
                    "cm_child_beta_p",
                    "cm_mother_beta",
                    "cm_mother_beta_se",
                    "cm_mother_beta_p",
                    "cf_cmf_p",
                    "cf_child_beta",
                    "cf_child_beta_se",
                    "cf_child_beta_p",
                    "cf_father_beta",
                    "cf_father_beta_se",
                    "cf_father_beta_p",
                    "mf_cmf_p",
                    "mf_mother_beta",
                    "mf_mother_beta_se",
                    "mf_mother_beta_p",
                    "mf_father_beta",
                    "mf_father_beta_se",
                    "mf_father_beta_p",
                    "c_cmf_p",
                    "c_cm_p",
                    "c_cf_p",
                    "c_beta",
                    "c_beta_se",
                    "c_beta_p",
                    "m_cmf_p",
                    "m_cm_p",
                    "m_mf_p",
                    "m_beta",
                    "m_beta_se",
                    "m_beta_p",
                    "f_cmf_p",
                    "f_cf_p",
                    "f_mf_p",
                    "f_beta",
                    "f_beta_se",
                    "f_beta_p"
            );
        String header = String.join(
                Utils.separator,
                "phenotype",
                "variantID",
                "h",
                "beta",
                "betaSE",
                "rSquare",
                "p",
                "nH",
                "n"
        );
        outputWriter.writeLine(header);

        try {

            ExecutorService pool = Executors.newFixedThreadPool(nVariants);

            IntStream.range(0, nVariants)
                    .mapToObj(
                            i -> new LinearModelRunnable(
                                    iterator,
                                    childToParentMap,
                                    phenotypesHandler,
                                    outputWriter
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
            iterator.close();

        }

        end = Instant.now().getEpochSecond();
        duration = end - start;

        System.out.println(
                Instant.now() + " - Done (" + iterator.getnVariants() + " variants processed in " + duration + " seconds)"
        );
    }
}
