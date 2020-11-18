package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.genome.ChromosomeUtils;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Generic iterator for a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfIteratorTargets implements VariantIterator {

    /**
     * The name of the file being iterated.
     */
    private final String fileName;
    /**
     * The variants to process.
     */
    private final VariantList variantList;
    /**
     * Variants with a distance to a target variant lower than maxDistance will
     * be included in the computation.
     */
    private final int maxDistance;
    /**
     * The file reader.
     */
    private final VCFFileReader reader;
    /**
     * Number of variant columns.
     */
    private int nVariantColumns;
    /**
     * Number of samples.
     */
    private int nSamples;
    /**
     * Mutex for the access to the file.
     */
    private final SimpleSemaphore mutex = new SimpleSemaphore(1);
    /**
     * The index in the variant list.
     */
    private int variantListIndex = -1;
    /**
     * The variant iterator.
     */
    private CloseableIterator<VariantContext> iterator = null;
    /**
     * The number of variants processed.
     */
    private int nVariants = 0;
    /**
     * The maximum number of variants to process, ignored if negative.
     */
    public static int nLimit = -1;
    /**
     * Boolean indicating whether the currently targeted variant has been found.
     */
    private boolean variantFound = false;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * Constructor.
     *
     * @param vcfFile The vcf file.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param variantList The ids of the variants to iterate.
     * @param logger A logger for issues processing variants.
     */
    public VcfIteratorTargets(
            File vcfFile,
            VariantList variantList,
            int maxDistance,
            SimpleCliLogger logger
    ) {

        File indexFile = getVcfIndexFile(vcfFile);
        fileName = vcfFile.getName();
        reader = new VCFFileReader(vcfFile, indexFile);
        this.variantList = variantList;
        this.maxDistance = maxDistance;
        this.logger = logger;

    }

    @Override
    public VcfVariant next() {

        mutex.acquire();

        if (iterator == null || !iterator.hasNext()) {

            if (iterator != null && variantListIndex < variantList.variantId.length) {

                System.out.println(
                        Instant.now() + " - " + variantList.variantId[variantListIndex] + " (" + (variantListIndex + 1) + " of " + variantList.variantId.length + ")"
                );

            }

            variantListIndex++;

            if (variantListIndex >= variantList.variantId.length) {

                mutex.release();
                return null;

            }
            
            if (variantListIndex >= 1 && !variantFound) {
                
                String variantId = variantList.variantId[variantListIndex - 1];
                        logger.logVariant(
                                variantId,
                                "Variant not found"
                        );
            }

            int windowStart = Math.max(variantList.start[variantListIndex] - maxDistance, 1);

            int targetBpEnd = variantList.end[variantListIndex] + maxDistance;
            Integer chrLength = ChromosomeUtils.chromosomeLength37.get(variantList.chromosome[variantListIndex]);

            int windowEnd = chrLength != null && chrLength < targetBpEnd ? chrLength : targetBpEnd;

            iterator = reader.query(
                    variantList.chromosome[variantListIndex],
                    windowStart,
                    windowEnd
            );
            
            variantFound = false;

        }

        VariantContext variantContext = iterator.next();

        if (variantContext == null) {

            mutex.release();
            return next();

        }

        String variantId = variantContext.getID();
        int variantBpStart = variantContext.getStart();
        int variantBpEnd = variantContext.getEnd();

        if (maxDistance == 0 && !variantId.equals(variantList.variantId[variantListIndex])
                || variantBpStart < variantList.start[variantListIndex] - maxDistance || variantBpEnd > variantList.end[variantListIndex] + maxDistance) {

            mutex.release();
            return next();

        } else if (variantId.equals(variantList.variantId[variantListIndex])) {
            
            variantFound = true;
            
        }

        nVariants++;

        mutex.release();

        return new VcfVariant(
                variantContext
        );

    }

    /**
     * Returns the number of variant columns.
     *
     * @return the number of variant columns
     */
    public int getnVariantColumns() {

        return nVariantColumns;

    }

    /**
     * Returns the number of samples.
     *
     * @return the number of samples
     */
    public int getnSamples() {

        return nSamples;

    }

    @Override
    public int getnVariants() {
        return nVariants;
    }

    @Override
    public void close() {

        iterator.close();;
        reader.close();

    }

    /**
     * Returns the index file for the given vcf file.
     *
     * @param vcfFile the vcf file
     *
     * @return the index file
     */
    public static File getVcfIndexFile(
            File vcfFile
    ) {

        return new File(vcfFile.getAbsolutePath() + ".tbi");

    }

    @Override
    public boolean isFinished() {

        return variantListIndex >= variantList.variantId.length;

    }

}
