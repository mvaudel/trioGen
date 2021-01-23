package no.uib.triogen.io.genotypes.bgen.writer;

import com.github.luben.zstd.Zstd;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.BitSet;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.BgenUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Writer for a bgen file with phased haplotypes.
 *
 * @author Marc Vaudel
 */
public class BgenWriter implements AutoCloseable {

    private final File destinationFile;
    /**
     * Random access file to the destination file.
     */
    private final RandomAccessFile raf;
    private final File indexFile;

    private final int headerBlockLength = 4 + 4 + 4 + 4 + 13 + 4;
    private int offset = -1;
    private int nVariants = -1;
    private String[] samples;

    private ArrayList<IndexEntry> indexEntries = new ArrayList<>();

    private final static int probabilityPrecision = 8;

    public BgenWriter(File destinationFile, File indexFile) throws IOException {

        this.destinationFile = destinationFile;
        this.indexFile = indexFile;

        this.raf = new RandomAccessFile(destinationFile, "rw");

    }

    public void addVariant(
            VariantInformation variantInformation,
            ArrayList<String[]> genotypedAlleles
    ) throws IOException {

        byte[] variantId = variantInformation.id.getBytes(ENCODING);
        raf.writeShort(Short.reverseBytes((short) variantId.length));
        raf.write(variantId);

        byte[] rsId = variantInformation.rsId.getBytes(ENCODING);
        raf.writeShort(Short.reverseBytes((short) rsId.length));
        raf.write(rsId);

        byte[] contig = variantInformation.contig.getBytes(ENCODING);
        raf.writeShort(Short.reverseBytes((short) contig.length));
        raf.write(contig);

        raf.writeInt(Integer.reverseBytes(variantInformation.position));

        String[] alleles = variantInformation.alleles;

        raf.writeShort(Short.reverseBytes((short) alleles.length));

        for (String allele : alleles) {

            byte[] alleleByte = allele.getBytes(ENCODING);
            raf.writeInt(Integer.reverseBytes(alleleByte.length));
            raf.write(alleleByte);

        }

        int minPloidy = -1;
        int maxPloidy = -1;
        int nProbabilities = 0;

        for (String[] genotypes : genotypedAlleles) {

            if (minPloidy == -1 || genotypes.length < minPloidy) {

                minPloidy = genotypes.length;

            }
            if (maxPloidy == -1 || genotypes.length > maxPloidy) {

                maxPloidy = genotypes.length;

            }

            nProbabilities += genotypes.length * (alleles.length - 1);

        }

        if (minPloidy < 0 || minPloidy > 63) {

            throw new IllegalArgumentException("Min ploidy of " + minPloidy + " for variant " + variantInformation.id + ", should be between 0 and 63 inclusive.");

        }
        if (maxPloidy < 0 || maxPloidy > 63) {

            throw new IllegalArgumentException("Max ploidy of " + maxPloidy + " for variant " + variantInformation.id + ", should be between 0 and 63 inclusive.");

        }

        byte[] probabilityData = new byte[4 + 2 + 1 + 1 + samples.length + 1 + 1 + nProbabilities]; // Note: only for precision of 8

        ByteBuffer buffer = ByteBuffer.wrap(probabilityData);

        buffer.putInt(Integer.reverseBytes(samples.length));
        buffer.putShort(Short.reverseBytes((short) alleles.length));
        buffer.put((byte) minPloidy);
        buffer.put((byte) maxPloidy);

        for (String[] genotypes : genotypedAlleles) {

            if (genotypes.length == 0) {

                buffer.put((byte) 0b10000000);

            } else {

                buffer.put((byte) genotypes.length);

            }
        }

        buffer.put((byte) 1);
        buffer.put((byte) probabilityPrecision);

        for (String[] genotypes : genotypedAlleles) {

            if (genotypes.length > 0) {

                for (int z = 0; z < genotypes.length; z++) {

                    for (int k = 0; k < alleles.length - 1; k++) {

                        if (alleles[k].equals(genotypes[z])) {

                            buffer.put((byte) 0b11111111);

                        } else {

                            buffer.put((byte) 0b00000000);

                        }
                    }
                }
            }
        }

        int maxLength = (int) Zstd.compressBound(probabilityData.length);

        byte[] destinationArray = new byte[maxLength];

        int compressedArrayLength = (int) Zstd.compress(destinationArray, probabilityData, 1);

        raf.writeInt(Integer.reverseBytes(compressedArrayLength + 4));

        long start = raf.getFilePointer();

        raf.writeInt(Integer.reverseBytes(probabilityData.length));

        raf.write(destinationArray, 0, compressedArrayLength);

        indexEntries.add(new IndexEntry(variantInformation, start, compressedArrayLength + 4));

    }

    public void initiate(String[] samples) throws IOException {

        this.samples = samples;

        raf.seek(headerBlockLength + 4);
        writeSampleIdentifierBlock();

        raf.seek(0);
        raf.writeInt(Integer.reverseBytes(offset));

        writeHeaderBlock();

        raf.seek(offset + 4);

    }

    public void finalize() throws IOException {

        raf.seek(8);
        raf.writeInt(Integer.reverseBytes(indexEntries.size()));

        try (SimpleFileWriter writer = new SimpleFileWriter(indexFile, true)) {

            writer.writeLine(BgenIndex.FIRST_LINE);

            writer.writeLine("# Compression:", Integer.toString(2));

            writer.writeLine("# Samples");

            writer.writeLine(samples);

            writer.writeLine("# Variants:", Integer.toString(indexEntries.size()));

            writer.writeLine("variantId", "rsId", "contig", "bp", "alleles", "index", "blockSize");

            for (int variantI = 0; variantI < indexEntries.size(); variantI++) {

                IndexEntry indexEntry = indexEntries.get(variantI);

                VariantInformation variantInformation = indexEntry.variantInformation;
                String id = variantInformation.id;
                long index = indexEntry.index;
                long length = indexEntry.blockSize;

                writer.writeLine(
                        id,
                        variantInformation.rsId,
                        variantInformation.contig,
                        Integer.toString(variantInformation.position),
                        String.join(",", variantInformation.alleles),
                        Long.toString(index),
                        Long.toString(length)
                );
            }
        }
    }

    private void writeHeaderBlock() throws IOException {

        raf.writeInt(Integer.reverseBytes(headerBlockLength));

        raf.writeInt(Integer.reverseBytes(nVariants));

        raf.writeInt(Integer.reverseBytes(samples.length));

        raf.write(BgenUtils.MAGIC_NUMBER);

        raf.write(BgenUtils.IDENTIFIER);

        BitSet bitSet = new BitSet(32);
        bitSet.set(1);
        bitSet.set(3);
        bitSet.set(31);

        raf.write(bitSet.toByteArray());

    }

    private void writeSampleIdentifierBlock() throws IOException {

        int nBytes = 8 + 2 * samples.length;

        ArrayList<byte[]> identifiersByte = new ArrayList<>(samples.length);

        for (String sample : samples) {

            byte[] identifierByte = sample.getBytes(ENCODING);

            identifiersByte.add(identifierByte);

            nBytes += identifierByte.length;

        }

        raf.writeInt(Integer.reverseBytes(nBytes));

        raf.writeInt(Integer.reverseBytes(samples.length));

        for (byte[] identifierByte : identifiersByte) {

            raf.writeShort(Short.reverseBytes((short) identifierByte.length));

            raf.write(identifierByte);

        }

        offset = (int) raf.getFilePointer() - 4;

    }

    @Override
    public void close() throws IOException {

        raf.close();

    }

    private class IndexEntry {

        public final VariantInformation variantInformation;

        public final long index;

        public final long blockSize;

        public IndexEntry(
                VariantInformation variantInformation,
                long index,
                long blockSize
        ) {

            this.variantInformation = variantInformation;
            this.index = index;
            this.blockSize = blockSize;

        }
    }
}
