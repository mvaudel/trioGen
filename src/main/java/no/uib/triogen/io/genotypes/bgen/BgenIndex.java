package no.uib.triogen.io.genotypes.bgen;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map.Entry;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Index for a bgen file.
 *
 * @author Marc Vaudel
 */
public class BgenIndex {

    private final static String FIRST_LINE = "# TrioGen_bgen_index_v.1.0";

    public final HashMap<String, VariantInformation> variantInformationMap;

    public final HashMap<String, Long> variantIndexMap;

    public final HashMap<String, Integer> sampleIndexMap;

    public BgenIndex(
            HashMap<String, VariantInformation> variantInformationMap,
            HashMap<String, Long> variantIndexMap,
            HashMap<String, Integer> sampleIndexMap
    ) {

        this.variantInformationMap = variantInformationMap;
        this.variantIndexMap = variantIndexMap;
        this.sampleIndexMap = sampleIndexMap;

    }

    public static BgenIndex getBgenIndex(
            File bgenFile
    ) throws IOException {

        return getBgenIndex(
                bgenFile,
                getDefaultIndexFile(bgenFile)
        );

    }

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

    public static void writeToFile(
            BgenIndex bgenIndex,
            File indexFile
    ) {

        try (SimpleFileWriter writer = new SimpleFileWriter(indexFile, true)) {

            writer.writeLine(FIRST_LINE);

            writer.writeLine("# Samples");

            String[] sampleIds = new String[bgenIndex.sampleIndexMap.size()];

            for (Entry<String, Integer> entry : bgenIndex.sampleIndexMap.entrySet()) {

                sampleIds[entry.getValue()] = entry.getKey();

            }

            writer.writeLine(sampleIds);

            writer.writeLine("# Variants:", Integer.toString(bgenIndex.variantInformationMap.size()));

            writer.writeLine("variantId", "rsId", "contig", "bp", "alleles", "index");

            for (VariantInformation variantInformation : bgenIndex.variantInformationMap.values()) {

                String id = variantInformation.id;

                long index = bgenIndex.variantIndexMap.get(id);

                writer.writeLine(
                        id,
                        variantInformation.rsId,
                        variantInformation.contig,
                        Integer.toString(variantInformation.bp),
                        String.join(",", variantInformation.alleles),
                        Long.toString(index)
                );
            }
        }
    }

    public static BgenIndex readFromFile(
            File indexFile
    ) {

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(indexFile)) {

            String line = reader.readLine();

            if (!line.equals(FIRST_LINE)) {

                throw new IllegalArgumentException("Index file " + indexFile + " could not be parsed as TrioGen index file. Please remove the file or provide a new location to save the index file.");

            }

            reader.readLine();
            line = reader.readLine();

            String[] samples = line.split(IoUtils.SEPARATOR);

            HashMap<String, Integer> sampleIndexMap = new HashMap<>(samples.length);

            for (int i = 0; i < samples.length; i++) {

                sampleIndexMap.put(samples[i], i);

            }

            line = reader.readLine();

            String[] lineSplit = line.split(IoUtils.SEPARATOR);
            int nVariants = Integer.parseInt(lineSplit[1]);

            HashMap<String, VariantInformation> variantInformationMap = new HashMap<>(nVariants);
            HashMap<String, Long> variantIndexMap = new HashMap<>(nVariants);

            reader.readLine();

            while ((line = reader.readLine()) != null) {

                if (line.length() > 0) {

                    lineSplit = line.split(IoUtils.SEPARATOR);

                    String id = lineSplit[0];
                    String rsId = lineSplit[1];
                    String contig = lineSplit[2];
                    int bp = Integer.parseInt(lineSplit[3]);
                    String[] alleles = lineSplit[4].split(",");
                    long position = Long.parseLong(lineSplit[5]);

                    VariantInformation variantInformation = new VariantInformation(id, rsId, contig, bp, alleles);

                    variantInformationMap.put(id, variantInformation);
                    variantIndexMap.put(id, position);

                }
            }

            return new BgenIndex(variantInformationMap, variantIndexMap, sampleIndexMap);

        }
    }

    public static File getDefaultIndexFile(
            File bgenFile
    ) {

        return new File(bgenFile.getAbsolutePath() + ".ti");

    }

    private static BgenIndex buidBgenIndex(
            File bgenFile
    ) throws IOException {

        HashMap<String, Integer> sampleIndexMap;
        HashMap<String, VariantInformation> variantInformationMap;
        HashMap<String, Long> variantIndexMap;

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

            if (!bitSet.get(0) || bitSet.get(1)) {

                throw new IllegalArgumentException("Only files with compressed SNP block probability are supported.");

            }

            if (bitSet.get(2) || !bitSet.get(3) || bitSet.get(4)) {

                throw new IllegalArgumentException("Only files with layout 2 are supported.");

            }

            if (!bitSet.get(31)) {

                throw new IllegalArgumentException("Only files with stored sample identifiers are supported.");

            }

            // Sample identifier
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

            sampleIndexMap = new HashMap<>(nSamples);

            for (int sampleI = 0; sampleI < nSamples; sampleI++) {

                short tempShort = Short.reverseBytes(raf.readShort());
                int sampleIdByteLength = Short.toUnsignedInt(tempShort);

                byte[] sampleIdBytes = new byte[sampleIdByteLength];
                raf.read(sampleIdBytes);

                String sampleId = new String(sampleIdBytes, 0, sampleIdByteLength, ENCODING);

                sampleIndexMap.put(sampleId, sampleI);

            }

            // Variant data block
            variantInformationMap = new HashMap<>(nVariants);
            variantIndexMap = new HashMap<>(nVariants);

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

                long bp = Integer.toUnsignedLong(tempInt);

                if (bp > Integer.MAX_VALUE || bp < 0) {

                    throw new IllegalArgumentException("Unexpected variant position (" + bp + ") for variant " + variantId + ", should be an integer between 0 (included) and " + Integer.MAX_VALUE + " (included).");

                }

                int bpInt = (int) bp;

                tempShort = Short.reverseBytes(raf.readShort());
                int nAlleles = Short.toUnsignedInt(tempShort);

                String[] alleles = new String[nAlleles];

                for (int alleleI = 0; alleleI < nAlleles; alleleI++) {

                    tempShort = Short.reverseBytes(raf.readShort());
                    int alleleLength = Short.toUnsignedInt(tempShort);

                    byte[] alleleBytes = new byte[alleleLength];
                    raf.read(contigBytes);

                    String allele = new String(alleleBytes, 0, alleleLength, ENCODING);

                    alleles[alleleI] = allele;

                }

                VariantInformation variantInformation = new VariantInformation(variantId, rsId, contig, bpInt, alleles);

                variantInformationMap.put(variantId, variantInformation);

                long index = raf.getFilePointer();

                variantIndexMap.put(variantId, index);

                tempInt = Integer.reverseBytes(raf.readInt());

                long variantDataLength = Integer.toUnsignedLong(tempInt);

                if (variantDataLength > Integer.MAX_VALUE || variantDataLength < 0) {

                    throw new IllegalArgumentException("Variant data length " + variantDataLength + " for variant " + variantId + " not supported, should be between 0 (included) and " + Integer.MAX_VALUE + " (included).");

                }

                toSkip = (int) variantDataLength;

                nSkipped = raf.skipBytes(toSkip);

                if (nSkipped != toSkip) {

                    throw new IllegalArgumentException("Unexpected number of bytes skipped in variant (expected: " + toSkip + "; skipped: " + nSkipped + ").");

                }
            }
        }

        return new BgenIndex(variantInformationMap, variantIndexMap, sampleIndexMap);

    }

}
