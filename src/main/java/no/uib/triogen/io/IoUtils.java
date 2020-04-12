package no.uib.triogen.io;

import java.io.File;
import java.lang.reflect.Method;
import java.nio.MappedByteBuffer;

/**
 * Utilities for reading and writing files.
 *
 * @author Marc Vaudel
 */
public class IoUtils {

    /**
     * Default encoding, cf the second rule.
     */
    public static final String ENCODING = "UTF-8";
    /**
     * Default separator
     */
    public static final String SEPARATOR = "\t";
    /**
     * The line separator.
     */
    public static final String LINE_SEPARATOR = System.getProperty("line.separator");

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
     * Attempts at closing a buffer to avoid memory issues. Adapted from
     * https://stackoverflow.com/questions/2972986/how-to-unmap-a-file-from-memory-mapped-using-filechannel-in-java.
     *
     * @param mappedByteBuffer The buffer to close.
     */
    public static void closeBuffer(
            MappedByteBuffer mappedByteBuffer
    ) {

        if (mappedByteBuffer == null || !mappedByteBuffer.isDirect()) {
            return;
        }

        try {

            Method cleaner = mappedByteBuffer.getClass().getMethod("cleaner");
            cleaner.setAccessible(true);
            Method clean = Class.forName("sun.misc.Cleaner").getMethod("clean");
            clean.setAccessible(true);
            clean.invoke(cleaner.invoke(mappedByteBuffer));

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
