package no.uib.triogen.utils;

import com.github.luben.zstd.Zstd;

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

        byte[] uncompressedByteAray = new byte[(int) uncompressedLength];

        long decompressedBytes = Zstd.decompress(uncompressedByteAray, compressedByteArray);

        if (decompressedBytes != uncompressedLength) {

            throw new IllegalArgumentException(decompressedBytes + " bytes decompressed where " + uncompressedLength + " expected.");

        }

        return uncompressedByteAray;

    }

    /**
     * Compresses the given byte array.
     *
     * @param uncompressedData The uncompressed data.
     *
     * @return The compressed data.
     */
    public static TempByteArray compress(
            byte[] uncompressedData
    ) {

        int maxLength = (int) Zstd.compressBound(uncompressedData.length);

        byte[] destinationArray = new byte[maxLength];

        int compressedArrayLength = (int) Zstd.compress(destinationArray, uncompressedData, 1);

        return new TempByteArray(destinationArray, compressedArrayLength);

    }

}
