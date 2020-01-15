package no.uib.triogen.io.flat.indexed;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import static no.uib.triogen.io.IoUtils.encoding;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class reads an indexed gz file. Note that unless otherwise specified, io
 * exceptions are thrown as runtime exceptions.
 *
 * @author Marc Vaudel
 */
public class IndexedGzReader implements AutoCloseable {

    /**
     * The random access file.
     */
    private final RandomAccessFile raf;
    /**
     * The inflater.
     */
    private final Inflater inflater = new Inflater(true);
    /**
     * Semaphore to synchronize threads.
     */
    private final SimpleSemaphore mutex = new SimpleSemaphore(1);

    /**
     * Constructor.
     * 
     * @param file The file to read from.
     * 
     * @throws IOException Exception thrown if an I/O error occurred.
     */
    public IndexedGzReader(
            File file
    ) throws IOException {

        raf = new RandomAccessFile(file, "r");

    }

    /**
     * Reads and uncompresses the content of the file at the given coordinates.
     * 
     * @param position The position in the file.
     * @param compressedLength The number of bytes to read.
     * @param uncompressedLength The number of bytes to uncompress.
     * 
     * @return the uncompressed content of the file.
     */
    public String read(
            long position,
            int compressedLength,
            int uncompressedLength
    ) {
        
        mutex.acquire();

        try {

            byte[] compressedByteArray = new byte[compressedLength];
            byte[] uncompressedByteAray = new byte[uncompressedLength];

            raf.seek(position);
            int bytesRead = raf.read(compressedByteArray);

            if (bytesRead != compressedLength) {

                throw new IllegalArgumentException("Unexpected number of bytes read " + bytesRead + " (expected: " + compressedLength + ")");

            }

            inflater.setInput(compressedByteArray);
            int bytesUncompressed = inflater.inflate(uncompressedByteAray);

            if (bytesUncompressed == 0) {

                throw new IllegalArgumentException("Missing input or dictionary.");

            } else if (bytesUncompressed != uncompressedLength) {

                throw new IllegalArgumentException("Unexpected number of bytes uncompressed " + bytesUncompressed + " (expected: " + uncompressedLength + ")");

            }

            return new String(uncompressedByteAray, 0, uncompressedByteAray.length, encoding);

        } catch (IOException e) {

            throw new RuntimeException(e);

        } catch (DataFormatException e) {

            throw new RuntimeException(e);

        } finally {
            
            mutex.release();
            
        }
    }

    @Override
    public void close() throws Exception {

        inflater.end();
        raf.close();

    }

}
