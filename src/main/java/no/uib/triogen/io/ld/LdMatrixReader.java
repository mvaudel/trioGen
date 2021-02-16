package no.uib.triogen.io.ld;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.mapping.MemoryMappedFile;
import no.uib.triogen.model.ld.R2;
import static no.uib.triogen.utils.CompressionUtils.zstdDecompress;

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
     * rsid to variant id map.
     */
    private final HashMap<String, ArrayList<String>> rsidToVariantIdMap;
    /**
     * variant id to rsid map.
     */
    private final HashMap<String, String> variantIdtoRsIdMap;
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

            byte[] footer = zstdDecompress(
                    compressedFooter, 
                    uncompressedLength
            );

            ByteBuffer byteBuffer = ByteBuffer.wrap(footer);
            int idsByteLength = byteBuffer.getInt();
            byte[] idsBytes = new byte[idsByteLength];
            byteBuffer.get(idsBytes);
            String idsString = new String(idsBytes, 0, idsByteLength, ENCODING);
            String[] allIds = idsString.split(IoUtils.SEPARATOR);
            int nIds = allIds.length / 2;
            variantIds = new String[nIds];
            variantIdtoRsIdMap = new HashMap<>(nIds);
            rsidToVariantIdMap = new HashMap<>(nIds);

            for (int i = 0; i < nIds; i++) {

                String id = allIds[i];
                String rsId = allIds[i + nIds];

                variantIds[i] = id;

                if (!rsId.equals("")) {

                    variantIdtoRsIdMap.put(id, rsId);

                    ArrayList<String> variantIdsForRsId = rsidToVariantIdMap.get(rsId);

                    if (variantIdsForRsId == null) {

                        variantIdsForRsId = new ArrayList<>(1);
                        rsidToVariantIdMap.put(rsId, variantIdsForRsId);

                    }

                    variantIdsForRsId.add(id);

                }
            }

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
     * @param variantA The id of variant A. It can be the id or the rsid.
     *
     * @return The id to r2 map.
     */
    public ArrayList<R2> getR2(
            String variantA
    ) {

        Long index = indexMap.get(variantA);

        if (index == null) {

            ArrayList<String> tempIds = getVariantIds(variantA); // See if it can be an rsid

            if (tempIds != null) {

                ArrayList<R2> result = new ArrayList<>(0);

                for (String variantId : tempIds) {

                    ArrayList<R2> tempResult = getR2(variantId);

                    if (tempResult != null) {

                        result.addAll(tempResult);

                    }
                }

                if (!result.isEmpty()) {

                    return result;

                }
            }

            return null;

        }

        int nVariants;
        byte[] compressedData;

        try (MemoryMappedFile.MiniBuffer buffer = memoryMappedFile.getBuffer(index)) {

            nVariants = buffer.getInt();

            if (nVariants == 0) {

                return new ArrayList<>(0);

            }

            int compressedLength = buffer.getInt();

            compressedData = new byte[compressedLength];
            buffer.get(compressedData);

        }

        int uncompressedLength = nVariants * Integer.BYTES + 2 * nVariants * Short.BYTES + nVariants * Float.BYTES;

        byte[] uncompressedData = zstdDecompress(compressedData, uncompressedLength);
        ByteBuffer byteBuffer = ByteBuffer.wrap(uncompressedData);

        ArrayList<R2> r2s = new ArrayList<>(nVariants);

        for (int i = 0; i < nVariants; i++) {

            int variantIndex = byteBuffer.getInt();
            short alleleA = byteBuffer.getShort();
            short alleleB = byteBuffer.getShort();
            float r2 = byteBuffer.getFloat();

            r2s.add(
                    new R2(
                            variantIndex,
                            alleleA,
                            alleleB,
                            r2
                    )
            );
        }

        return r2s;

    }

    /**
     * Returns the id corresponding to a given variant index.
     *
     * @param variantIndex The index of the variant.
     *
     * @return The id of the variant.
     */
    public String getId(
            int variantIndex
    ) {

        return variantIds[variantIndex];

    }

    /**
     * Returns the rsid for the given variant, null if none set.
     *
     * @param variantId The id of the variant.
     *
     * @return The rsid of the variant.
     */
    public String getRsId(
            String variantId
    ) {

        return variantIdtoRsIdMap.get(variantId);

    }

    /**
     * Returns the variant ids for the given rsId, null if none set.
     *
     * @param rsid The rsid of the variant.
     *
     * @return The variant ids.
     */
    public ArrayList<String> getVariantIds(
            String rsid
    ) {

        return rsidToVariantIdMap.get(rsid);

    }

    @Override
    public void close() throws Exception {

        memoryMappedFile.close();

    }
}
