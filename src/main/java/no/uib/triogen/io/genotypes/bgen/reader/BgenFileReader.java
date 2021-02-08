package no.uib.triogen.io.genotypes.bgen.reader;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * Reader for a Bgen file.
 *
 * @author Marc Vaudel
 */
public class BgenFileReader implements AutoCloseable {

    private final BgenIndex bgenIndex;
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
    private final MappedByteBuffer[] mappedByteBuffers;
    /**
     * The index of the variant in the bgen file.
     */
    private final int[] variantIndexes;

    private final HashMap<Integer, char[]> inheritanceMap;

    public BgenFileReader(
            File bgenFile,
            BgenIndex bgenIndex,
            VariantList variantList,
            int distance,
            HashMap<Integer, char[]> inheritanceMap
    ) throws IOException {

        this.bgenIndex = bgenIndex;
        this.inheritanceMap = inheritanceMap;

        raf = new RandomAccessFile(bgenFile, "r");

        fc = raf.getChannel();

        if (variantList == null) {

            variantIndexes = null;
            mappedByteBuffers = new MappedByteBuffer[bgenIndex.variantIdArray.length];

            for (int i = 0; i < bgenIndex.variantIdArray.length; i++) {

                long blockStart = bgenIndex.variantIndexArray[i];
                long blockLength = bgenIndex.variantBlockLengthArray[i];

                mappedByteBuffers[i] = fc.map(
                        FileChannel.MapMode.READ_ONLY,
                        blockStart,
                        blockLength
                );

            }

        } else {

            ArrayList<MappedByteBuffer> bufferList = new ArrayList<>(variantList.variantId.length);
            ArrayList<Integer> indexesList = new ArrayList<>(variantList.variantId.length);

            for (int i = 0; i < bgenIndex.variantIdArray.length; i++) {

                VariantInformation variantInformation = bgenIndex.variantInformationArray[i];

                boolean found = false;

                for (int j = 0; j < variantList.variantId.length; j++) {

                    if (variantInformation.position >= variantList.position[j] - distance
                            && variantInformation.position <= variantList.position[j] + distance
                            && variantInformation.contig.equals(variantList.contig[j])) {

                        found = true;
                        break;

                    }
                }

                if (found) {

                    long blockStart = bgenIndex.variantIndexArray[i];
                    long blockLength = bgenIndex.variantBlockLengthArray[i];

                    MappedByteBuffer mappedByteBuffer = fc.map(
                            FileChannel.MapMode.READ_ONLY,
                            blockStart,
                            blockLength
                    );

                    bufferList.add(mappedByteBuffer);
                    indexesList.add(i);

                }
            }

            mappedByteBuffers = bufferList.toArray(new MappedByteBuffer[bufferList.size()]);

            variantIndexes = indexesList.stream()
                    .mapToInt(
                            i -> i
                    )
                    .toArray();

        }
    }

    public BgenFileReader(
            File bgenFile,
            BgenIndex bgenIndex,
            HashMap<Integer, char[]> inheritanceMap
    ) throws IOException {

        this(bgenFile, bgenIndex, null, 0, inheritanceMap);
        
    }

    public int getNVariants() {

        return mappedByteBuffers.length;

    }

    public VariantInformation getVariantInformation(int i) {

        int index = variantIndexes == null ? i : variantIndexes[i];

        return bgenIndex.variantInformationArray[index];

    }

    public ByteBuffer getDataBlock(int i) {

        int index = variantIndexes == null ? i : variantIndexes[i];

        MappedByteBuffer mappedByteBuffer = mappedByteBuffers[index];
        return mappedByteBuffer.duplicate();

    }

    public long getBlockLength(int i) {

        int index = variantIndexes == null ? i : variantIndexes[i];

        return bgenIndex.variantBlockLengthArray[index];

    }

    public int nVariants() {

        return mappedByteBuffers.length;

    }

    public BgenVariantData getVariantData(int i) {

        VariantInformation variantInformation = getVariantInformation(i);
        ByteBuffer buffer = getDataBlock(i);
        long blockLength = getBlockLength(i);

        if (blockLength > Integer.MAX_VALUE) {

            throw new IllegalArgumentException("Block length (" + blockLength + ") for variant " + variantInformation.id + " exceeds maximal capacity (" + Integer.MAX_VALUE + ").");

        }

        return new BgenVariantData(
                bgenIndex.sampleIds,
                variantInformation,
                buffer,
                (int) blockLength,
                bgenIndex.compressionType,
                inheritanceMap
        );
    }

    @Override
    public void close() throws Exception {

        raf.close();
        fc.close();

        for (MappedByteBuffer mappedByteBuffer : mappedByteBuffers) {

            IoUtils.closeBuffer(mappedByteBuffer);

        }

    }
}
