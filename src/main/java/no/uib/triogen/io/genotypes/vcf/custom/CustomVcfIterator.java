package no.uib.triogen.io.genotypes.vcf.custom;

import java.io.File;
import java.time.Instant;
import java.util.HashMap;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This iterator iterates through a vcf file formatted according to the output of the Sanger imputation server.
 *
 * @author Marc Vaudel
 */
public class CustomVcfIterator implements VariantIterator {

    /**
     * The name of the file being iterated.
     */
    private final String fileName;
    /**
     * The file reader.
     */
    private final SimpleFileReader reader;
    /**
     * Sample to column index map.
     */
    private final HashMap<String, Integer> sampleMap;
    /**
     * Number of variant columns.
     */
    private int nVariantColumns;
    /**
     * Number of samples.
     */
    private int nSamples;
    /**
     * Boolean indicating whether the end of file was reached.
     */
    private boolean endOfFile = false;
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
     * Boolean indicating whether the genotypes provider should cache genotype
     * values.
     */
    private final boolean useCache;

    /**
     * Constructor.
     *
     * @param file The file to iterate.
     * @param useCache Boolean indicating whether the genotypes provider should
     * cache genotype values.
     */
    public CustomVcfIterator(
            File file,
            boolean useCache
    ) {
        
        this.useCache = useCache;

        // Set up reader
        fileName = file.getName();
        reader = SimpleFileReader.getFileReader(file);

        // Move to header
        String line;
        while ((line = reader.readLine()) != null && line.charAt(1) == '#') {

        }
        if (line == null) {

            throw new IllegalArgumentException("Reached end of file when looking for the header.");

        }

        // Parse header
        String[] lineSplit = line.split("\t");

        nVariantColumns = lineSplit[8].equals("FORMAT") ? 9 : 8;
        nSamples = lineSplit.length - nVariantColumns;
        sampleMap = new HashMap<>(nSamples);

        for (int i = nVariantColumns; i < lineSplit.length; i++) {

            sampleMap.put(lineSplit[i], i - nVariantColumns);

        }
    }

    @Override
    public VcfLine next() {

        mutex.acquire();

        String line = endOfFile ? null : reader.readLine();
        nVariants++;
        endOfFile = line == null || nLimit != -1 && nVariants >= nLimit;

        mutex.release();

        if (!endOfFile) {

            if (nVariants >= progress + 100000) {

                progress = nVariants;

                System.out.println(
                        Instant.now() + " - " + fileName + " " + nVariants + " variants processed"
                );

            }
        }

        return line == null ? null
                : new VcfLine(
                        this,
                        line,
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

    /**
     * Returns the column index for a sample.
     *
     * @param sampleId the id of the sample
     *
     * @return the column index for a sample
     */
    public int getSampleIndex(
            String sampleId
    ) {

        return sampleMap.get(sampleId);

    }

    @Override
    public int getnVariants() {
        return nVariants;
    }

    @Override
    public void close() {

        reader.close();

    }

    @Override
    public boolean isFinished() {
        
        return endOfFile;
        
    }
}
