package no.uib.triogen.io.flat.indexed.gz;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.zip.CRC32;
import java.util.zip.Deflater;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.utils.SimpleSemaphore;
import uk.ac.ebi.pride.tools.braf.BufferedRandomAccessFile;

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
    public final int headerLength = 10;
    /**
     * GZIP header magic number.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private final static int GZIP_MAGIC = 0x8b1f;

    /**
     * Compression method for the deflate algorithm (the only one currently
     * supported).
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    public static final int DEFLATED = 8;
    /**
     * CRC-32 of uncompressed data.
     *
     * Adapted from java.util.zip.GZIPOutputStream by David Connelly. Copyright
     * (c) 1996, 2013, Oracle and/or its affiliates. No copyright infringement
     * intended.
     */
    private CRC32 crc = new CRC32();
    /**
     * The deflater used to compress data.
     */
    private final Deflater deflater;
    /**
     * The random access file used to write data.
     */
    private final RandomAccessFile bw;
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

        bw = new RandomAccessFile(file, "rw");

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

        bw.write(new byte[]{
            (byte) GZIP_MAGIC, // Magic number (short)
            (byte) (GZIP_MAGIC >> 8), // Magic number (short)
            DEFLATED, // Compression method (CM)
            0, // Flags (FLG)
            0, // Modification time MTIME (int)
            0, // Modification time MTIME (int)
            0, // Modification time MTIME (int)
            0, // Modification time MTIME (int)
            0, // Extra flags (XFLG)
            0 // Operating system (OS)
        });
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

            byte[] inputBytes = inputString.getBytes(IoUtils.encoding);
            byte[] output = new byte[inputBytes.length];

            deflater.setInput(inputBytes);

            int compressedDataLength = deflater.deflate(output, 0, output.length, Deflater.FULL_FLUSH);

            if (compressedDataLength > 0) {

                bw.write(output, 0, compressedDataLength);

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

            byte[] trailer = new byte[8];

            writeInt((int) crc.getValue(), trailer, 0); // CRC-32 of uncompr. data
            writeInt(deflater.getTotalIn(), trailer, 4); // Number of uncompr. bytes
            bw.write(trailer);

            deflater.end();
            bw.close();

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
}
