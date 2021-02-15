package no.uib.triogen.io.genotypes.bgen.reader;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
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

    /**
     * The index of the bgen file.
     */
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
     * The allele inheritance map.
     */
    private final HashMap<Integer, char[]> inheritanceMap;
    /**
     * The default ploidy for mothers.
     */
    private final int defaultMotherPloidy;
    /**
     * The default ploidy for fathers.
     */
    private final int defaultFatherPloidy;

    /**
     * Constructor.
     *
     * @param bgenFile The file to read.
     * @param bgenIndex The index of the file.
     * @param inheritanceMap The allele inheritance map to use.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     *
     * @throws IOException Exception thrown if an error occurs while reading the
     * file.
     */
    public BgenFileReader(
            File bgenFile,
            BgenIndex bgenIndex,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy
    ) throws IOException {

        this.bgenIndex = bgenIndex;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;

        raf = new RandomAccessFile(bgenFile, "r");

        fc = raf.getChannel();

    }

    /**
     * Returns information on the given variant.
     *
     * @param i The index of the variant of interest.
     *
     * @return Information on the given variant.
     */
    public VariantInformation getVariantInformation(int i) {

        return bgenIndex.variantInformationArray[i];

    }

    /**
     * Returns a buffer wrapped around the data block of the given variant.
     *
     * @param i The index of the variant of interest.
     *
     * @return A buffer wrapped around the data block of the given variant.
     */
    public MappedByteBuffer getDataBlock(int i) {

        try {

                long blockStart = bgenIndex.variantIndexArray[i];
                long blockLength = bgenIndex.variantBlockLengthArray[i];

                return fc.map(
                        FileChannel.MapMode.READ_ONLY,
                        blockStart,
                        blockLength
                );

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Returns the length of the block of the given variant.
     *
     * @param i The index of the variant of interest.
     *
     * @return The length of the block of the given variant.
     */
    public long getBlockLength(int i) {

        return bgenIndex.variantBlockLengthArray[i];

    }

    /**
     * Returns the variant data for the given variant.
     *
     * @param i The index of the variant of interest.
     *
     * @return The variant data for the given variant.
     */
    public BgenVariantData getVariantData(int i) {

        VariantInformation variantInformation = getVariantInformation(i);
        MappedByteBuffer buffer = getDataBlock(i);
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
                inheritanceMap,
                defaultMotherPloidy,
                defaultFatherPloidy
        );
    }

    @Override
    public void close() throws Exception {

        raf.close();
        fc.close();
        
    }
}
