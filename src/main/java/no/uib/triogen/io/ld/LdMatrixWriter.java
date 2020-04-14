package no.uib.triogen.io.ld;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.zip.Deflater;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.utils.TempByteArray;
import static no.uib.triogen.io.ld.LdMatrixUtils.MAGIC_NUMBER;
import no.uib.triogen.model.geno.VariantIndex;
import static no.uib.triogen.utils.CompressionUtils.compress;
import no.uib.triogen.utils.SimpleSemaphore;
import no.uib.triogen.utils.Utils;

/**
 * Writer for an ld matrix.
 *
 * @author Marc Vaudel
 */
public class LdMatrixWriter implements AutoCloseable {

    /**
     * The length of the file header.
     */
    public static final int HEADER_LENGTH = MAGIC_NUMBER.length + Long.BYTES;
    /**
     * The random access file to write to.
     */
    private final RandomAccessFile raf;
    /**
     * Index for the variants.
     */
    private final VariantIndex variantIndex;
    /**
     * List of the variant indexes.
     */
    private final ArrayList<Integer> variantIndexes = new ArrayList<>();
    /**
     * List of the index in the file.
     */
    private final ArrayList<Long> indexesInFile = new ArrayList<>();
    /**
     * Semaphore to synchronize threads writing to the file.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param outputFile the output file.
     * @param variantIndex The index to use for the variants.
     *
     * @throws FileNotFoundException Exception thrown if the output file was not
     * found.
     * @throws IOException Exception thrown if an error occurred while
     * attempting to write to output file.
     */
    public LdMatrixWriter(
            VariantIndex variantIndex,
            File outputFile
    ) throws FileNotFoundException, IOException {

        if (outputFile.exists()) {

            outputFile.delete();

        }

        this.variantIndex = variantIndex;

        raf = new RandomAccessFile(outputFile, "rw");
        raf.seek(HEADER_LENGTH);

    }

    /**
     * Writes a variant to the file.
     *
     * @param variantIndex The index of the variant.
     * @param variantIds The indexes of the other variants.
     * @param r2s The ld r2s between the variant and the other variants.
     * @param deflater The deflater to compress parts of the file.
     *
     * @throws IOException Exception thrown if an error occurred while
     * attempting to write to output file.
     */
    public void addVariant(
            int variantIndex,
            ArrayList<Integer> variantIds,
            ArrayList<Double> r2s,
            Deflater deflater
    ) throws IOException {

        int nVariants = variantIds.size();
        ByteBuffer buffer = ByteBuffer.allocate(nVariants * Integer.BYTES + nVariants * Double.BYTES);

        for (int i = 0; i < nVariants; i++) {

            buffer.putInt(variantIds.get(i));
            buffer.putDouble(r2s.get(i));

        }

        byte[] uncompressedData = buffer.array();
        TempByteArray compressedData = compress(uncompressedData, deflater);
        
        if (compressedData.length == 0) {
            
            throw new IllegalArgumentException("Empty compressed data.");
            
        }

        semaphore.acquire();

        long index = raf.getFilePointer() - HEADER_LENGTH;

        variantIndexes.add(variantIndex);
        indexesInFile.add(index);

        raf.writeInt(nVariants);
        raf.writeInt(compressedData.length);
        raf.write(compressedData.array, 0, compressedData.length);

        semaphore.release();

    }

    /**
     * Compresses and writes the given byte array.
     *
     * @param uncompressedData The uncompressed data.
     * @param deflater The deflater to compress parts of the file.
     *
     * @throws IOException Exception thrown if an error occurred while
     * attempting to write to output file.
     */
    private void compressAndWrite(
            byte[] uncompressedData,
            Deflater deflater
    ) throws IOException {

        TempByteArray compressedData = compress(uncompressedData, deflater);

        raf.writeInt(compressedData.length);
        raf.writeInt(uncompressedData.length);
        raf.write(compressedData.array, 0, compressedData.length);

    }

    /**
     * Writes the header and footer to the file.
     *
     * @throws IOException Exception thrown if an error occurred while
     * attempting to write to output file.
     */
    private void writeHeaderAndFooter() throws IOException {

        long footerPosition = raf.getFilePointer();

        String[] variantIds = variantIndex.getVariants();

        String variantIdsString = Arrays.stream(variantIds)
                .collect(
                        Collectors.joining(SEPARATOR)
                );

        byte[] titleBytes = variantIdsString.getBytes(IoUtils.ENCODING);

        ByteBuffer buffer = ByteBuffer.allocate(2 * Integer.BYTES + titleBytes.length + indexesInFile.size() * Integer.BYTES + indexesInFile.size() * Long.BYTES);

        buffer
                .putInt(titleBytes.length)
                .put(titleBytes)
                .putInt(variantIndexes.size());

        for (int i = 0; i < variantIndexes.size(); i++) {

            int variantI = variantIndexes.get(i);
            buffer.putInt(variantI);

            long index = indexesInFile.get(i);
            buffer.putLong(index);

        }

        Deflater deflater = new Deflater(Deflater.BEST_COMPRESSION, true);
        compressAndWrite(buffer.array(), deflater);

        raf.seek(0);
        raf.write(MAGIC_NUMBER);
        raf.writeLong(footerPosition);

    }

    @Override
    public void close() throws IOException {

        writeHeaderAndFooter();

        raf.close();

    }

}
