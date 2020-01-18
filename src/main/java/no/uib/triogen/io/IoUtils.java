package no.uib.triogen.io;

import java.io.File;
import java.lang.reflect.Method;
import java.nio.MappedByteBuffer;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.readers.SimpleGzReader;
import no.uib.triogen.io.flat.readers.SimpleTextReader;

/**
 * Utilities for reading and writing files.
 *
 * @author Marc Vaudel
 */
public class IoUtils {

    /**
     * Default encoding, cf the second rule.
     */
    public static final String encoding = "UTF-8";
    /**
     * Default separator
     */
    public static final String separator = "\t";
    /**
     * The line separator.
     */
    public static final String lineSeparator = System.getProperty("line.separator");

    /**
     * Returns the index file for the given vcf file.
     *
     * @param vcfFilePath the complete vcf file path
     *
     * @return the index file
     */
    public static File getVcfIndexFile(File vcfFilePath) {

        return new File(vcfFilePath + ".tbi");

    }

    /**
     * Returns a file reader for the given file. Gz reader if the file ends with
     * ".gz", text reader otherwise.
     *
     * @param file the file to read
     *
     * @return a file reader for the given file
     */
    public static SimpleFileReader getFileReader(File file) {

        return file.getName().endsWith(".gz") ? new SimpleGzReader(file) : new SimpleTextReader(file);

    }

    /**
     * Attempts at closing a buffer to avoid memory issues. Adapted from
     * https://stackoverflow.com/questions/2972986/how-to-unmap-a-file-from-memory-mapped-using-filechannel-in-java.
     *
     * @param buffer The buffer to close.
     */
    public static void closeBuffer(
            MappedByteBuffer buffer
    ) {

        if (buffer == null || !buffer.isDirect()) {
            return;
        }

        try {

            Method cleaner = buffer.getClass().getMethod("cleaner");
            cleaner.setAccessible(true);
            Method clean = Class.forName("sun.misc.Cleaner").getMethod("clean");
            clean.setAccessible(true);
            clean.invoke(cleaner.invoke(buffer));

        } catch (Exception ex) {
        }
    }

    /**
     * Returns the index corresponding to the given file.
     *
     * @param file The file indexed.
     *
     * @return The index file.
     */
    public static File getIndexFile(
            File file
    ) {

        String stem = file.getAbsolutePath();

        if (stem.endsWith(".gz")) {

            stem = stem.substring(0, stem.length() - 3);

        }

        return new File(stem + ".index.gz");
    }
}
