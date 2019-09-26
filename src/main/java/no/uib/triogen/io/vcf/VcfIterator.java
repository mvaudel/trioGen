package no.uib.triogen.io.vcf;

import java.io.File;
import java.util.HashMap;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 *
 * @author Marc Vaudel
 */
public class VcfIterator implements AutoCloseable {

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
     * The maximum number of variants to process, ignored if negative.
     */
    public static int nLimit = -1;

    /**
     * Constructor.
     *
     * @param file the file.
     */
    public VcfIterator(File file) {

        // Set up reader
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

    /**
     * Reads the next line as unparsed VCF line. Returns null if reading the
     * file is finished.
     *
     * @return the next line
     */
    public VcfLine next() {

        mutex.acquire();

        String line = endOfFile ? null : reader.readLine();
        endOfFile = line == null || nLimit != -1 && ++nVariants > nLimit;

        mutex.release();

        return line == null ? null
                : new VcfLine(
                        this,
                        line
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
    public int getSampleIndex(String sampleId) {

        return sampleMap.get(sampleId);

    }

    /**
     * Returns the number of variants read from the file.
     * 
     * @return the number of variants read from the file
     */
    public int getnVariants() {
        return nVariants;
    }
    
    

    @Override
    public void close() {

        reader.close();

    }
}
