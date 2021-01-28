package no.uib.cell_rk.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import no.uib.triogen.utils.SimpleSemaphore;
import static no.uib.triogen.io.IoUtils.ENCODING;
import static no.uib.triogen.io.IoUtils.SEPARATOR;

/**
 * Simple gz file writer that is thread safe and throws exceptions as runtime
 * exceptions.
 *
 * @author Marc Vaudel
 */
public class SimpleFileWriter implements AutoCloseable {

    /**
     * The buffered writer used for writing.
     */
    private final BufferedWriter bw;
    /**
     * A simple mutex.
     */
    private final SimpleSemaphore mutex = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param file the file to write to
     * @param gz boolean indicating whether the output should be gzipped
     */
    public SimpleFileWriter(
            File file, 
            boolean gz
    ) {

        try {

            if (gz) {

                FileOutputStream fileStream = new FileOutputStream(file);
                GZIPOutputStream gzipStream = new GZIPOutputStream(fileStream);
                OutputStreamWriter encoder = new OutputStreamWriter(gzipStream, ENCODING);
                bw = new BufferedWriter(encoder);

            } else {

                bw = new BufferedWriter(new FileWriter(file));

            }

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Writes a new line using the give elements. Elements are separated using
     * the separator in the properties class.
     *
     * @param elements line elements
     */
    public void writeLine(
            String... elements
    ) {

        String line = String.join(SEPARATOR, elements);
        writeLine(line);

    }

    /**
     * Writes a new line.
     *
     * @param line the line to write
     */
    public void writeLine(
            String line
    ) {

        write(line, true);

    }

    public void newLine() {

        try {

            bw.newLine();

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Writes some text.
     *
     * @param text the text to write
     * @param newLine boolean indicating whether an end of line should be
     * appended
     */
    public void write(
            String text, 
            boolean newLine
    ) {

        try {

            mutex.acquire();

            bw.write(text);

            if (newLine) {
                bw.newLine();
            }

            mutex.release();

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    @Override
    public void close() {

        try {

            bw.close();

        } catch (IOException e) {

            throw new RuntimeException(e);

        }

    }

}
