package no.uib.triogen.utils;

import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Functions needed for compression and decompression.
 *
 * @author Marc Vaudel
 */
public class CompressionUtils {
    

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
     * Compresses the given byte array.
     *
     * @param uncompressedData The uncompressed data.
     * @param deflater The deflater to compress parts of the file.
     *
     * @return The compressed data.
     */
    public static TempByteArray compress(
            byte[] uncompressedData,
            Deflater deflater
    ) {

        byte[] compressedData = new byte[uncompressedData.length];

        int outputLength = compressedData.length;

        deflater.setInput(uncompressedData);
        int compressedByteLength = deflater.deflate(
                compressedData,
                0,
                compressedData.length,
                Deflater.FULL_FLUSH
        );
        int compressedDataLength = compressedByteLength;

        while (compressedByteLength == outputLength) {

            byte[] output2 = new byte[outputLength];
            compressedByteLength = deflater.deflate(
                    output2,
                    0,
                    outputLength,
                    Deflater.FULL_FLUSH
            );

            compressedData = Utils.mergeArrays(
                    compressedData,
                    output2,
                    compressedByteLength
            );
            compressedDataLength += compressedByteLength;

        }

        return new TempByteArray(compressedData, compressedDataLength);

    }
    
    
}
