package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
import no.uib.triogen.model.geno.VariantList;
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
     * Boolean indicating whether the genotypes provider should cache genotype
     * values.
     */
    private final boolean useCache;

    /**
     * Constructor.
     *
     * @param vcfFile The vcf file-
     * @param variantList The ids of the variants to iterate.
     * @param useCache Boolean indicating whether the genotypes provider should
     * cache genotype values.
     */
    public VcfIteratorTargets(
            File vcfFile,
            VariantList variantList,
            boolean useCache
    ) {

        this.useCache = useCache;

        File indexFile = getVcfIndexFile(vcfFile);
        fileName = vcfFile.getName();
        reader = new VCFFileReader(vcfFile, indexFile);
        this.variantList = variantList;

    }

    @Override
    public VcfVariant next() {

        mutex.acquire();

        if (iterator == null || !iterator.hasNext()) {

            variantListIndex++;

            if (variantListIndex >= variantList.variantId.length) {

                mutex.release();
                return null;

            }

            iterator = reader.query(
                    variantList.chromosome[variantListIndex],
                    variantList.start[variantListIndex],
                    variantList.end[variantListIndex]
            );

        }

        VariantContext variantContext = iterator.next();

        if (variantContext == null) {

            mutex.release();
            return next();

        }

        String vcfVariantId = variantContext.getID();

        if (vcfVariantId == null || !vcfVariantId.equals(variantList.variantId[variantListIndex])) {

            mutex.release();
            return next();

        }

        System.out.println(
                Instant.now() + " - " + variantList.variantId[variantListIndex] + " (" + (variantListIndex + 1) + " of " + variantList.variantId.length + ")"
        );

        nVariants++;

        mutex.release();

        return new VcfVariant(
                variantContext,
                useCache
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

}
