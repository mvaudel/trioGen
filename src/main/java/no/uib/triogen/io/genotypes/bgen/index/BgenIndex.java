package no.uib.triogen.io.genotypes.bgen.index;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.BitSet;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.BgenUtils;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Index for a bgen file.
 *
 * @author Marc Vaudel
 */
public class BgenIndex {

    /**
     * The first line of the index.
     */
    public final static String FIRST_LINE = "# TrioGen_bgen_index_v.1.0.1";
/**
 * Array of the variant ids.
 */
    public final String[] variantIdArray;
/**
 * Array of the information on the variants.
 */
    public final VariantInformation[] variantInformationArray;
/**
 * Array of the indexes.
 */
    public final long[] variantIndexArray;
/**
 * Array of the lengths of the blocks.
 */
    public final long[] variantBlockLengthArray;
/**
 * Array of the ids of the samples.
 */
    public final String[] sampleIds;
    /**
     * The compression type according to the format specifications.
     */
    public final int compressionType;

    /**
     * Constructor.
     * 
     * @param variantIdArray Array of the variant ids.
     * @param variantInformationArray Array of the information on the variants.
     * @param variantIndexArray Array of the indexes.
     * @param variantBlockLengthArray Array of the lengths of the blocks.
     * @param sampleIds Array of the ids of the samples.
     * @param compressionType The compression type according to the format specifications.
     */
    public BgenIndex(
            String[] variantIdArray,
            VariantInformation[] variantInformationArray,
            long[] variantIndexArray,
            long[] variantBlockLengthArray,
            String[] sampleIds,
            int compressionType
    ) {

        this.variantIdArray = variantIdArray;
        this.variantInformationArray = variantInformationArray;
        this.variantIndexArray = variantIndexArray;
        this.variantBlockLengthArray = variantBlockLengthArray;
        this.sampleIds = sampleIds;
        this.compressionType = compressionType;

    }

    /**
     * Returns an index for the given bgen file.
     * 
     * @param bgenFile The bgen file to index.
     * 
     * @return The index.
     * 
     * @throws IOException Exception thrown if an error occurred while reading or writing a file.
     */
    public static BgenIndex getBgenIndex(
            File bgenFile
    ) throws IOException {

        return getBgenIndex(
                bgenFile,
                getDefaultIndexFile(bgenFile)
        );

    }

    /**
     * Returns the bgen index for the given file.
     * 
     * @param bgenFile The bgen file to index.
     * @param indexFile The index file where the index should be saved.
     * 
     * @return The index.
     * 
     * @throws IOException Exception thrown if an error occurred while reading or writing a file.
     */
    public static BgenIndex getBgenIndex(
            File bgenFile,
            File indexFile
    ) throws IOException {

        if (indexFile.exists()) {

            return readFromFile(indexFile);

        } else {

            BgenIndex bgenIndex = buidBgenIndex(bgenFile);

            writeToFile(bgenIndex, indexFile);

            return bgenIndex;

        }
    }

    /**
     * Writes the given index to the given file.
     * 
     * @param bgenIndex The bgen index.
     * @param indexFile The file where to save the index.
     */
    public static void writeToFile(
            BgenIndex bgenIndex,
            File indexFile
    ) {

        try (SimpleFileWriter writer = new SimpleFileWriter(indexFile, true)) {

            writer.writeLine(FIRST_LINE);

            writer.writeLine("# Compression:", Integer.toString(bgenIndex.compressionType));

            writer.writeLine("# Samples");

            writer.writeLine(bgenIndex.sampleIds);

            writer.writeLine("# Variants:", Integer.toString(bgenIndex.variantIdArray.length));

            writer.writeLine("variantId", "rsId", "contig", "bp", "alleles", "index", "blockSize");

            for (int variantI = 0; variantI < bgenIndex.variantIdArray.length; variantI++) {

                String id = bgenIndex.variantIdArray[variantI];
                VariantInformation variantInformation = bgenIndex.variantInformationArray[variantI];
                long start = bgenIndex.variantIndexArray[variantI];
                long length = bgenIndex.variantBlockLengthArray[variantI];

                writer.writeLine(id,
                        variantInformation.rsid,
                        variantInformation.contig,
                        Integer.toString(variantInformation.position),
                        String.join(",", variantInformation.alleles),
                        Long.toString(start),
                        Long.toString(length)
                );
            }
        }
    }

    /**
     * Reads the bgen index from the given file.
     * 
     * @param indexFile The file where the index is saved.
     * 
     * @return The index.
     */
    public static BgenIndex readFromFile(
            File indexFile
    ) {

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(indexFile, false)) {

            String line = reader.readLine();

            if (!line.equals(FIRST_LINE)) {

                throw new IllegalArgumentException("Index file " + indexFile + " could not be parsed as TrioGen index file. Please remove the file or provide a new location to save the index file.");

            }

            line = reader.readLine();
            String[] lineSplit = line.split(IoUtils.SEPARATOR);
            int compressionType = Integer.parseInt(lineSplit[1]);

            reader.readLine();
            line = reader.readLine();
            String[] samples = line.split(IoUtils.SEPARATOR);

            line = reader.readLine();
            lineSplit = line.split(IoUtils.SEPARATOR);
            int nVariants = Integer.parseInt(lineSplit[1]);

            String[] variantIdArray = new String[nVariants];
            VariantInformation[] variantInformationArray = new VariantInformation[nVariants];
            long[] variantIndexArray = new long[nVariants];
            long[] variantBlockLengthArray = new long[nVariants];

            reader.readLine();

            int variantI = 0;

            while ((line = reader.readLine()) != null) {

                if (line.length() > 0) {

                    lineSplit = line.split(IoUtils.SEPARATOR);

                    String id = lineSplit[0];
                    String rsId = lineSplit[1];
                    String contig = lineSplit[2];
                    int bp = Integer.parseInt(lineSplit[3]);
                    String[] alleles = lineSplit[4].split(",");
                    long start = Long.parseLong(lineSplit[5]);
                    long blockSize = Long.parseLong(lineSplit[6]);

                    VariantInformation variantInformation = new VariantInformation(id, rsId, contig, bp, alleles);

                    variantInformationArray[variantI] = variantInformation;
                    variantIndexArray[variantI] = start;
                    variantBlockLengthArray[variantI] = blockSize;
                    variantIdArray[variantI] = id;

                    variantI++;

                }
            }

            return new BgenIndex(variantIdArray, variantInformationArray, variantIndexArray, variantBlockLengthArray, samples, compressionType);

        }
    }

    /**
     * Returns the default index file for the given bgen file.
     * 
     * @param bgenFile The bgen file.
     * 
     * @return The default file where to save the index.
     */
    public static File getDefaultIndexFile(
            File bgenFile
    ) {

        return new File(bgenFile.getAbsolutePath() + ".index.gz");

    }

    /**
     * Builds the bgen index for the given bgen file.
     * 
     * @param bgenFile The bgen file.
     * 
     * @return The bgen index.
     * 
     * @throws IOException Exception thrown if an error occurred while reading or writing a file.
     */
    private static BgenIndex buidBgenIndex(
            File bgenFile
    ) throws IOException {

        String[] samples;
        String[] variantIdArray;
        VariantInformation[] variantInformationArray;
        long[] variantIndexArray;
        long[] variantBlockLengthArray;
        int compressionType;

        try (RandomAccessFile raf = new RandomAccessFile(bgenFile, "r")) {

            // Offset
            int tempInt = Integer.reverseBytes(raf.readInt());

            long offSet = Integer.toUnsignedLong(tempInt);

            // Header
            tempInt = Integer.reverseBytes(raf.readInt());

            long headerBlockLength = Integer.toUnsignedLong(tempInt);

            if (headerBlockLength > offSet) {

                throw new IllegalArgumentException("Header block (" + headerBlockLength + ") larger than offset (" + offSet + ") in bgen file " + bgenFile + ".");

            }
            if (headerBlockLength - 20 > Integer.MAX_VALUE) {

                throw new IllegalArgumentException("Header block exceeds maximal size (" + Integer.MAX_VALUE + " + 20).");

            }
            if (headerBlockLength - 20 < 0) {

                throw new IllegalArgumentException("Header block too small (should contain at least 20 bytes).");

            }

            tempInt = Integer.reverseBytes(raf.readInt());

            long nVariantsLong = Integer.toUnsignedLong(tempInt);

            if (nVariantsLong <= 0 || nVariantsLong >= Integer.MAX_VALUE) {

                throw new IllegalArgumentException("Unexpected number of variants (" + nVariantsLong + "), should be an integer between 0 (excluded) and " + Integer.MAX_VALUE + " (excluded).");

            }

            int nVariants = (int) nVariantsLong;

            tempInt = Integer.reverseBytes(raf.readInt());

            long nSamples1 = Integer.toUnsignedLong(tempInt);

            byte[] fileMagicNumber = new byte[4];
            raf.read(fileMagicNumber);

            if (!BgenUtils.checkMagicNumber(fileMagicNumber)) {

                throw new IllegalArgumentException("Magic number (" + new String(fileMagicNumber) + ") should be 'bgen' or '0000'.");

            }

            long otherDataLength = (int) headerBlockLength - 20;

            if (otherDataLength > Integer.MAX_VALUE || otherDataLength < 0) {

                throw new IllegalArgumentException("Header data length " + otherDataLength + " not supported, should be between 0 (included) and " + Integer.MAX_VALUE + " (included).");

            }

            int toSkip = (int) otherDataLength;

            int nSkipped = raf.skipBytes(toSkip);

            if (nSkipped != toSkip) {

                throw new IllegalArgumentException("Unexpected number of bytes skipped in header (expected: " + toSkip + "; skipped: " + nSkipped + ").");

            }

            byte[] flags = new byte[4];
            raf.read(flags);

            BitSet bitSet = BitSet.valueOf(flags);
            
            if (!bitSet.get(0) && bitSet.get(1)) {
                
                compressionType = 2;
                
            } else if (bitSet.get(0) && !bitSet.get(1)) {
                
                compressionType = 1;
                
            } else if (!bitSet.get(0) && !bitSet.get(1)) {
                
                compressionType = 0;
                
            } else {
                
                throw new IllegalArgumentException("Compression value " + bitSet + " not supported.");
                
            }

            if (bitSet.get(2) || !bitSet.get(3) || bitSet.get(4)) {

                throw new IllegalArgumentException("Only files with layout 2 are supported.");

            }

            if (!bitSet.get(31)) {

                throw new IllegalArgumentException("Only files with stored sample identifiers are supported.");

            }

            // Sample identifiers
            tempInt = Integer.reverseBytes(raf.readInt());

            long identifierBlockLength = Integer.toUnsignedLong(tempInt);

            if (identifierBlockLength + headerBlockLength > offSet) {

                throw new IllegalArgumentException("Identifier block length larger than expected (identifier block length: " + identifierBlockLength + "; header block length: " + headerBlockLength + "; offset: " + offSet + ").");

            }

            tempInt = Integer.reverseBytes(raf.readInt());

            long nSamples2 = Integer.toUnsignedLong(tempInt);

            if (nSamples1 != nSamples2) {

                throw new IllegalArgumentException("Number of samples mismatch between header and identifier blocks (header: " + nSamples1 + "; identifier: " + nSamples2 + ").");

            }

            if (nSamples1 <= 0 || nSamples1 > Integer.MAX_VALUE) {

                throw new IllegalArgumentException("Unexpected number of samples (" + nSamples1 + "), should be an integer between 0 (excluded) and " + Integer.MAX_VALUE + " (included).");

            }

            int nSamples = (int) nSamples1;
            
            samples = new String[nSamples];

            for (int sampleI = 0; sampleI < nSamples; sampleI++) {

                short tempShort = Short.reverseBytes(raf.readShort());
                int sampleIdByteLength = Short.toUnsignedInt(tempShort);

                byte[] sampleIdBytes = new byte[sampleIdByteLength];
                raf.read(sampleIdBytes);

                String sampleId = new String(sampleIdBytes, 0, sampleIdByteLength, ENCODING);
                
                samples[sampleI] = sampleId;

            }

            // Variant data block
            variantIdArray = new String[nVariants];
            variantInformationArray = new VariantInformation[nVariants];
            variantIndexArray = new long[nVariants];
            variantBlockLengthArray = new long[nVariants];

            for (int variantI = 0; variantI < nVariants; variantI++) {

                // Variant identifying data
                short tempShort = Short.reverseBytes(raf.readShort());
                int variantIdByteLength = Short.toUnsignedInt(tempShort);

                byte[] variantIdBytes = new byte[variantIdByteLength];
                raf.read(variantIdBytes);

                String variantId = new String(variantIdBytes, 0, variantIdByteLength, ENCODING);

                tempShort = Short.reverseBytes(raf.readShort());
                int rsIdByteLength = Short.toUnsignedInt(tempShort);

                byte[] rsIdBytes = new byte[rsIdByteLength];
                raf.read(rsIdBytes);

                String rsId = new String(rsIdBytes, 0, rsIdByteLength, ENCODING);

                tempShort = Short.reverseBytes(raf.readShort());
                int contigByteLength = Short.toUnsignedInt(tempShort);

                byte[] contigBytes = new byte[contigByteLength];
                raf.read(contigBytes);

                String contig = new String(contigBytes, 0, contigByteLength, ENCODING);

                tempInt = Integer.reverseBytes(raf.readInt());

                long bpLong = Integer.toUnsignedLong(tempInt);

                if (bpLong > Integer.MAX_VALUE || bpLong < 0) {

                    throw new IllegalArgumentException("Unexpected variant position (" + bpLong + ") for variant " + variantId + ", should be an integer between 0 (included) and " + Integer.MAX_VALUE + " (included).");

                }

                int bp = (int) bpLong;

                tempShort = Short.reverseBytes(raf.readShort());
                int nAlleles = Short.toUnsignedInt(tempShort);

                String[] alleles = new String[nAlleles];

                for (int alleleI = 0; alleleI < nAlleles; alleleI++) {

                    tempInt = Integer.reverseBytes(raf.readInt());

                    long alleleLengthLong = Integer.toUnsignedLong(tempInt);

                    if (alleleLengthLong > Integer.MAX_VALUE || alleleLengthLong < 0) {

                        throw new IllegalArgumentException("Allele length " + alleleLengthLong + " for allele " + alleleI + " of variant " + variantId + " not supported, should be between 0 (included) and " + Integer.MAX_VALUE + " (included).");

                    }

                    int alleleLength = (int) alleleLengthLong;

                    byte[] alleleBytes = new byte[alleleLength];
                    raf.read(alleleBytes);

                    String allele = new String(alleleBytes, 0, alleleLength, ENCODING);

                    alleles[alleleI] = allele;

                }

                if (variantId.length() == 0) {

                    variantId = String.join("_",
                            contig,
                            Integer.toString(bp),
                            String.join("_",
                                    alleles
                            )
                    );
                }

                variantIdArray[variantI] = variantId;

                VariantInformation variantInformation = new VariantInformation(variantId, rsId, contig, bp, alleles);

                variantInformationArray[variantI] = variantInformation;

                tempInt = Integer.reverseBytes(raf.readInt());

                long blockSizeLong = Integer.toUnsignedLong(tempInt);

                if (blockSizeLong > Integer.MAX_VALUE || blockSizeLong < 0) {

                    throw new IllegalArgumentException("Block size " + blockSizeLong + " for variant " + variantId + " not supported, should be between 0 (included) and " + Integer.MAX_VALUE + " (included).");

                }

                int blockSize = (int) blockSizeLong;

                long index = raf.getFilePointer();

                variantIndexArray[variantI] = index;
                variantBlockLengthArray[variantI] = blockSize;

                toSkip = blockSize;

                nSkipped = raf.skipBytes(toSkip);

                if (nSkipped != toSkip) {

                    throw new IllegalArgumentException("Unexpected number of bytes skipped in variant (expected: " + toSkip + "; skipped: " + nSkipped + ").");

                }
            }
        }

        return new BgenIndex(variantIdArray, variantInformationArray, variantIndexArray, variantBlockLengthArray, samples, compressionType);

    }
}
