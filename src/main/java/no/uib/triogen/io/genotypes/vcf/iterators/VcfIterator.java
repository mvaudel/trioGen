package no.uib.triogen.io.genotypes.vcf.iterators;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.Closeable;
import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.vcf.reader.VcfVariant;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Generic iterator for a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfIterator implements Closeable {

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
     * The sample ids.
     */
    private String[] sampleIds;
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
        
        sampleIds = reader.getFileHeader().getSampleNamesInOrder().stream().toArray(String[]::new);

    }

    /**
     * Returns the next variant in the vcf file, null if none.
     * 
     * @return The next variant in the vcf file, null if none.
     */
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
     * Returns the sample ids.
     *
     * @return the sample ids
     */
    public String[] getSamplesIds() {

        return sampleIds;

    }

    /**
     * Returns the number of variants processed.
     * 
     * @return The number of variants processed.
     */
    public int getnVariants() {
        return nVariants;
    }

    @Override
    public void close() {

        iterator.close();
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

    /**
     * Returns a boolean indicating whether the iteration is finished, ie there is no next variant.
     * 
     * @return A boolean indicating whether the iteration is finished, ie there is no next variant.
     */
    public boolean isFinished() {
        
        return !iterator.hasNext();
        
    }

}
