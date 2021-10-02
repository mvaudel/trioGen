package no.uib.triogen.utils;

import io.airlift.compress.zstd.ZstdCompressor;
import io.airlift.compress.zstd.ZstdDecompressor;

/**
 * Functions needed for compression and decompression.
 *
 * @author Marc Vaudel
 */
public class CompressionUtils {

    /**
     * Decompresses the given byte array.
     *
     * @param compressedByteArray The compressed byte array.
     * @param uncompressedLength The uncompressed length.
     *
     * @return The decompressed array.
     */
    public static byte[] zstdDecompress(
            byte[] compressedByteArray,
            int uncompressedLength
    ) {

        return zstdDecompress(
                new ZstdDecompressor(),
                compressedByteArray,
                uncompressedLength
        );
    }

    /**
     * Decompresses the given byte array.
     *
     * @param decompressor The decompressor to use.
     * @param compressedByteArray The compressed byte array.
     * @param uncompressedLength The uncompressed length.
     *
     * @return The decompressed array.
     */
    public static byte[] zstdDecompress(
            ZstdDecompressor decompressor,
            byte[] compressedByteArray,
            int uncompressedLength
    ) {

        byte[] uncompressedByteAray = new byte[(int) uncompressedLength];

        long decompressedBytes = decompressor.decompress(compressedByteArray, 0, compressedByteArray.length, uncompressedByteAray, 0, uncompressedLength);

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
    public static TempByteArray zstdCompress(
            byte[] uncompressedData
    ) {

        return zstdCompress(
                new ZstdCompressor(),
                uncompressedData
        );

    }

    /**
     * Compresses the given byte array.
     *
     * @param compressor The compressor to use.
     * @param uncompressedData The uncompressed data.
     *
     * @return The compressed data.
     */
    public static TempByteArray zstdCompress(
            ZstdCompressor compressor,
            byte[] uncompressedData
    ) {

        if (uncompressedData.length == 0) {

            return new TempByteArray(new byte[0], 0);

        }

        int maxLength = (int) compressor.maxCompressedLength(uncompressedData.length);

        if (maxLength <= 0) {

            maxLength = uncompressedData.length;

        }

        byte[] destinationArray = new byte[maxLength];

        int compressedArrayLength = (int) compressor.compress(
                uncompressedData,
                0,
                uncompressedData.length,
                destinationArray,
                0,
                maxLength
        );

        return new TempByteArray(destinationArray, compressedArrayLength);

    }

}
