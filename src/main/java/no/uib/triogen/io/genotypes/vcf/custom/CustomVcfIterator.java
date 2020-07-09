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
     * The sample ids as they appear in the vcf file.
     */
    public final String[] samples;
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
     * Constructor.
     *
     * @param file The file to iterate.
     */
    public CustomVcfIterator(
            File file
    ) {

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
        
        samples = new String[nSamples];
        
        System.arraycopy(
                lineSplit, 
                nVariantColumns, 
                samples, 
                0, 
                nSamples
        );
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
