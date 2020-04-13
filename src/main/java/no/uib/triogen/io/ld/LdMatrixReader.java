package no.uib.triogen.io.ld;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.mapping.MemoryMappedFile;
import static no.uib.triogen.utils.CompressionUtils.uncompress;

/**
 * Reader for an ld matrix.
 *
 * @author Marc Vaudel
 */
public class LdMatrixReader implements AutoCloseable {

    /**
     * The block size to use for the memory mapped file.
     */
    public static final long BLOCK_SIZE = 64 * 1024 * 1024;
    /**
     * The ids of the variants.
     */
    public final String[] variantIds;
    /**
     * The index of the variants.
     */
    private final HashMap<String, Long> indexMap;
    /**
     * The memory mapped file.
     */
    private final MemoryMappedFile memoryMappedFile;

    /**
     * Constructor.
     *
     * @param file The file to read.
     *
     * @throws FileNotFoundException Exception thrown if the file was not found.
     * @throws IOException Exception thrown if an error occurred while
     * attempting to read the file.
     */
    public LdMatrixReader(
            File file
    ) throws FileNotFoundException, IOException {

        RandomAccessFile raf = new RandomAccessFile(file, "r");

        try {

            byte[] fileMagicNumber = new byte[LdMatrixUtils.MAGIC_NUMBER.length];
            raf.read(fileMagicNumber);

            if (!Arrays.equals(LdMatrixUtils.MAGIC_NUMBER, fileMagicNumber)) {

                raf.close();
                throw new IOException("File format of " + file + " not supported.");

            }

            long footerPosition = raf.readLong();

            raf.seek(footerPosition);
            int compressedLength = raf.readInt();
            int uncompressedLength = raf.readInt();

            byte[] compressedFooter = new byte[compressedLength];
            raf.read(compressedFooter);

            byte[] footer = uncompress(compressedFooter, uncompressedLength);

            ByteBuffer byteBuffer = ByteBuffer.wrap(footer);

            int idsByteLength = byteBuffer.getInt();
            byte[] idsBytes = new byte[idsByteLength];
            byteBuffer.get(idsBytes);
            String idsString = new String(idsBytes, 0, idsByteLength, ENCODING);
            variantIds = idsString.split(IoUtils.SEPARATOR);

            int nVariants = byteBuffer.getInt();
            indexMap = new HashMap<>(nVariants);

            for (int i = 0; i < nVariants; i++) {

                int variantI = byteBuffer.getInt();
                String variantId = variantIds[variantI];
                long index = byteBuffer.getLong();
                indexMap.put(variantId, index);

            }

            long offset = LdMatrixWriter.HEADER_LENGTH;
            long length = footerPosition - offset;

            this.memoryMappedFile = new MemoryMappedFile(
                    file,
                    offset,
                    length,
                    BLOCK_SIZE
            );

        } finally {
            raf.close();
        }
    }

    /**
     * Returns a map containing the ids of the variants B in LD with variant A
     * and the corresponding r2.
     *
     * @param variantA The id of variant A.
     *
     * @return The id to r2 map.
     */
    public HashMap<String, Double> getR2(
            String variantA
    ) {

        Long index = indexMap.get(variantA);

        if (index == null) {

            return null;

        }
        
        if (index > BLOCK_SIZE) {
            int debug = 1;
        }

        int nVariants;
        byte[] compressedData;

        try ( MemoryMappedFile.MiniBuffer buffer = memoryMappedFile.getBuffer(index)) {

            nVariants = buffer.getInt();
            int compressedLength = buffer.getInt();
            
            compressedData = new byte[compressedLength];
            buffer.get(compressedData);

        }

        int uncompressedLength = nVariants * Integer.BYTES + nVariants * Double.BYTES;

        byte[] uncompressedData = uncompress(compressedData, uncompressedLength);
        ByteBuffer byteBuffer = ByteBuffer.wrap(uncompressedData);

        HashMap<String, Double> result = new HashMap<>(nVariants);

        for (int i = 0; i < nVariants; i++) {

            int variantIndex = byteBuffer.getInt();
            String variantB = variantIds[variantIndex];
            double r2 = byteBuffer.getDouble();

            result.put(variantB, r2);

        }

        return result;

    }

    @Override
    public void close() throws Exception {
        
        memoryMappedFile.close();
        
    }
}
