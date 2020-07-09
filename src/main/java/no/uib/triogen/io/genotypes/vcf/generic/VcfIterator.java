package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Generic iterator for a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfIterator implements VariantIterator {

    /**
     * The name of the file being iterated.
     */
    private final String fileName;
    /**
     * The file reader.
     */
    private final VCFFileReader reader;
    /**
     * The variant iterator.
     */
    private final CloseableIterator<VariantContext> iterator;
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
     * The number of variants processed.
     */
    private int nVariants = 0;
    /**
     * The last progress.
     */
    private int progress = 0;
    /**
     * The maximum number of variants to process, ignored if negative.
     */
    public static int nLimit = -1;

    /**
     * Constructor.
     *
     * @param vcfFile The vcf file.
     */
    public VcfIterator(
            File vcfFile
    ) {

        File indexFile = getVcfIndexFile(vcfFile);
        fileName = vcfFile.getName();
        reader = new VCFFileReader(vcfFile, indexFile);
        iterator = reader.iterator();

    }

    @Override
    public VcfVariant next() {

        mutex.acquire();

        if (!iterator.hasNext()) {

            mutex.release();
            return null;

        }

        VariantContext variantContext = iterator.next();

        nVariants++;

        mutex.release();

        if (nVariants >= progress + 100000) {

            progress = nVariants;

            System.out.println(
                    Instant.now() + " - " + fileName + " " + nVariants + " variants processed"
            );

        }

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
        
        return !iterator.hasNext();
        
    }

}
