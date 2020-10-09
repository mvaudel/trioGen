package no.uib.triogen.io.flat.indexed;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.zip.CRC32;
import java.util.zip.Deflater;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.utils.SimpleSemaphore;
import static no.uib.triogen.utils.Utils.mergeArrays;

/**
 * This class writes an indexed gz file. Note that unless otherwise specified,
 * io exceptions are thrown as runtime exceptions.
 *
 * Parts of this class are adapted from java.util.zip.GZIPOutputStream by David
 * Connelly. Copyright (c) 1996, 2013, Oracle and/or its affiliates. No
 * copyright infringement intended.
 *
 * @author Marc Vaudel
 */
public class IndexedGzWriter implements AutoCloseable {
    
    /**
     * Length of the header in bytes. Content starts at this position.
     */
    public static final int HEADER_LENGTH = 10;
    /**
     * GZIP header magic number.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private final static int GZIP_MAGIC = 0x8b1f;
    /**
     * CRC-32 of uncompressed data.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private final CRC32 crc = new CRC32();
    /**
     * The deflater used to compress data.
     */
    private final Deflater deflater;
    /**
     * The file being written.
     */
    private final File file;
    /**
     * The random access file used to write data.
     */
    private final RandomAccessFile raf;
    /**
     * A mutex to synchronize the writing.
     */
    private final SimpleSemaphore mutex = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param file The file to write to.
     * @param compressionLevel The compression level to use.
     *
     * @throws IOException Exception thrown if an i/o error occurs.
     */
    public IndexedGzWriter(
            File file,
            int compressionLevel
    ) throws IOException {
        
        this.file = file;

        raf = new RandomAccessFile(file, "rw");

        writeHeader();
        crc.reset();

        deflater = new Deflater(compressionLevel, true);

    }

    /**
     * Constructor.
     *
     * @param file The file to write to.
     *
     * @throws IOException Exception thrown if an i/o error occurs.
     */
    public IndexedGzWriter(
            File file
    ) throws IOException {

        this(file, Deflater.DEFAULT_COMPRESSION);

    }

    /**
     * Writes GZIP member header.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private void writeHeader() throws IOException {

        raf.write(
                new byte[]{
                    (byte) GZIP_MAGIC, // Magic number (short)
                    (byte) (GZIP_MAGIC >> 8), // Magic number (short)
                    Deflater.DEFLATED, // Compression method (CM)
                    0, // Flags (FLG)
                    0, // Modification time MTIME (int)
                    0, // Modification time MTIME (int)
                    0, // Modification time MTIME (int)
                    0, // Modification time MTIME (int)
                    0, // Extra flags (XFLG)
                    0 // Operating system (OS)
                }
        );
    }

    /**
     * Appends a String to the file and returns the coordinates of the file.
     * Note that the deflater is flushed after this input, so the performance of
     * the compression depends on the length of the input: shorter input will
     * yield decreased compression rate.
     *
     * @param inputString The input as String
     *
     * @return The coordinates of the input in the file.
     */
    public IndexedGzCoordinates append(String inputString) {

        mutex.acquire();

        try {

            byte[] inputBytes = inputString.getBytes(IoUtils.ENCODING);
            byte[] output = new byte[inputBytes.length];

            deflater.setInput(inputBytes);

            int outputLength = output.length;
            int compressedDataLength = deflater.deflate(output, 0, outputLength, Deflater.FULL_FLUSH);
            int compressedByteLength = compressedDataLength;

            while (compressedByteLength == outputLength) {

                byte[] output2 = new byte[outputLength];
                compressedByteLength = deflater.deflate(output2, 0, outputLength, Deflater.FULL_FLUSH);

                output = mergeArrays(output, output2, compressedByteLength);
                compressedDataLength += compressedByteLength;

            }

            if (compressedDataLength > 0) {

                raf.write(output, 0, compressedDataLength);

            }

            crc.update(inputBytes);

            return new IndexedGzCoordinates(
                    compressedDataLength,
                    inputBytes.length
            );

        } catch (IOException e) {

            throw new RuntimeException(e);

        } finally {

            mutex.release();

        }
    }

    @Override
    public void close() {

        try {

            deflater.finish();

            while (!deflater.finished()) {

                byte[] output = new byte[512];

                int outputLength = output.length;
                int compressedDataLength = deflater.deflate(output, 0, outputLength, Deflater.FULL_FLUSH);
                int compressedByteLength = compressedDataLength;

                while (compressedByteLength == outputLength) {

                    byte[] output2 = new byte[outputLength];
                    compressedByteLength = deflater.deflate(output2, 0, outputLength, Deflater.FULL_FLUSH);

                    output = mergeArrays(output, output2, compressedByteLength);
                    compressedDataLength += compressedByteLength;

                }

                if (compressedDataLength > 0) {

                    raf.write(output, 0, compressedDataLength);

                }
            }

            int crcValue = (int) crc.getValue();
            int deflaterInput = deflater.getTotalIn();
            deflater.end();

            byte[] trailer = new byte[8];
            writeInt(crcValue, trailer, 0); // CRC-32 of uncompr. data
            writeInt(deflaterInput, trailer, 4); // Number of uncompr. bytes

            raf.write(trailer);
            raf.close();

        } catch (Exception e) {

            throw new RuntimeException(e);

        }

    }

    /*
     * Writes integer in Intel byte order to a byte array, starting at a
     * given offset.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private void writeInt(
            int i,
            byte[] buf,
            int offset
    ) {

        writeShort(i & 0xffff, buf, offset);
        writeShort((i >> 16) & 0xffff, buf, offset + 2);

    }

    /*
     * Writes short integer in Intel byte order to a byte array, starting
     * at a given offset
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private void writeShort(
            int s,
            byte[] buf,
            int offset
    ) {

        buf[offset] = (byte) (s & 0xff);
        buf[offset + 1] = (byte) ((s >> 8) & 0xff);

    }

    /**
     * Returns the file being written.
     * 
     * @return The file being written.
     */
    public File getFile() {
        return file;
    }
}
