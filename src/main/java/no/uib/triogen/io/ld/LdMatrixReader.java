package no.uib.triogen.io.ld;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Reader for an ld matrix.
 *
 * @author Marc Vaudel
 */
public class LdMatrixReader {

    /**
     * The ids of the variants.
     */
    public final String[] variantIds;
    /**
     * The index of the variants.
     */
    private final HashMap<String, Integer> indexMap;
    /**
     * Semaphore to synchronize threads.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);
    /**
     * The random access file.
     */
    private final RandomAccessFile raf;
    /**
     * The channel to the file.
     */
    private final FileChannel fc;
    /**
     * The mapped byte buffer.
     */
    private final MappedByteBuffer mappedByteBuffer;
    
    /**
     * Constructor.
     * 
     * @param file The file to read.
     * 
     * @throws FileNotFoundException Exception thrown if the file was not found.
     * @throws IOException Exception thrown if an error occurred while attempting to read
     * the file.
     */
    public LdMatrixReader(
            File file
    ) throws FileNotFoundException, IOException {
        
        raf = new RandomAccessFile(file, "r");

        try {
            
            byte[] fileMagicNumber = new byte[LdMatrixUtils.MAGIC_NUMBER.length];
            raf.read(fileMagicNumber);

            if (!Arrays.equals(LdMatrixUtils.MAGIC_NUMBER, fileMagicNumber)) {

                raf.close();
                throw new IOException("File format of " + file + " not supported.");

            }

            long footerPosition = raf.readLong();

            raf.seek(footerPosition);
            int length = raf.readInt();
            int uncompressedLength = raf.readInt();

            byte[] compressedFooter = new byte[length];
            raf.read(compressedFooter);

            byte[] footer = uncompress(compressedFooter, uncompressedLength);
            
            ByteBuffer byteBuffer = ByteBuffer.wrap(footer);
            
            int idsByteLength = byteBuffer.getInt();
            byte[] idsBytes = new byte[idsByteLength];
            byteBuffer.get(idsBytes);
            String idsString = new String(idsBytes, 0, idsByteLength, ENCODING);
            variantIds = idsString.split(IoUtils.SEPARATOR);
            
            indexMap = new HashMap<>(variantIds.length);

            for (int i = 0; i < variantIds.length; i++) {

                int index = byteBuffer.getInt();
                indexMap.put(variantIds[i], index);

            }

            fc = raf.getChannel();

            long size = footerPosition - LdMatrixWriter.HEADER_LENGTH;

            mappedByteBuffer = fc.map(FileChannel.MapMode.READ_ONLY, LdMatrixWriter.HEADER_LENGTH, size);

        } finally {
            raf.close();
        }
    }

    /**
     * Uncompresses the given byte array.
     *
     * @param compressedByteArray the compressed byte array
     * @param uncompressedLength the uncompressed length
     *
     * @return the uncompressed array
     */
    public static byte[] uncompress(
            byte[] compressedByteArray,
            int uncompressedLength
    ) {

        try {

            byte[] uncompressedByteAray = new byte[uncompressedLength];

            Inflater inflater = new Inflater(true);

            inflater.setInput(compressedByteArray);
            int bytesUncompressed = inflater.inflate(uncompressedByteAray);

            if (bytesUncompressed == 0) {

                throw new IllegalArgumentException("Missing input or dictionary.");

            } else if (bytesUncompressed != uncompressedLength) {

//                String debug = new String(uncompressedByteAray, 0, uncompressedByteAray.length, encoding);
                throw new IllegalArgumentException("Unexpected number of bytes uncompressed " + bytesUncompressed + " (expected: " + uncompressedLength + ")");

            }

            return uncompressedByteAray;

        } catch (DataFormatException e) {

            throw new RuntimeException(e);

        }
    }
    
    /**
     * Returns a map containing the ids of the variants B in LD with variant A and the corresponding r2.
     * 
     * @param variantA The id of variant A.
     * 
     * @return The id to r2 map.
     */
    public HashMap<String, Double> getR2(
    String variantA
    ) {
        
        Integer index = indexMap.get(variantA);
        
        if (index == null) {
            
            return null;
            
        }
        
        semaphore.acquire();

        mappedByteBuffer.position(index);
        
        int nVariants = mappedByteBuffer.getInt();
        int compressedLength = mappedByteBuffer.getInt();
        byte[] compressedData = new byte[compressedLength];
        mappedByteBuffer.get(compressedData);

        semaphore.release();
        
        int uncompressedLength = 2 * nVariants * Integer.BYTES;

        byte[] uncompressedData = uncompress(compressedData, uncompressedLength);
        ByteBuffer byteBuffer = ByteBuffer.wrap(uncompressedData);
        
        HashMap<String, Double> result = new HashMap<>(nVariants);
        
        for (int i = 0 ; i < nVariants ; i++) {
            
            int variantIndex = byteBuffer.getInt();
            String variantB = variantIds[variantIndex];
            double r2 = byteBuffer.getDouble();
            
            result.put(variantB, r2);
            
        }
        
        return result;
        
    }
    
}
