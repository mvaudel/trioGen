package no.uib.triogen.io.genotypes.bgen.reader;

import com.github.luben.zstd.Zstd;
import java.nio.ByteBuffer;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.Inflater;
import static no.uib.triogen.io.genotypes.InheritanceUtils.FATHER;
import static no.uib.triogen.io.genotypes.InheritanceUtils.MOTHER;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Genotypes provider based on a bgen file.
 *
 * @author Marc Vaudel
 */
public class BgenVariantData {

    private final String[] sampleIds;
    private final VariantInformation variantInformation;
    private final ByteBuffer compressedDataBlockContent;
    private final int blockLength;
    /**
     * Values for the haplotype probabilities for the children and parents.
     */
    private HashMap<String, double[]> haplotypeProbabilities = null;
    private HashSet<String> missing = null;
    private int nAlleles = -1;
    private final int compressionType;
    private double[] alleleFrequencyChildren;
    private int nChildrenGenotyped;

    private final HashMap<Integer, char[]> inheritanceMap;

    public BgenVariantData(
            String[] sampleIds,
            VariantInformation variantInformation,
            ByteBuffer dataBlockContent,
            int blockLength,
            int compressionType,
            HashMap<Integer, char[]> inheritanceMap
    ) {

        this.sampleIds = sampleIds;
        this.variantInformation = variantInformation;
        this.compressedDataBlockContent = dataBlockContent;
        this.blockLength = blockLength;
        this.compressionType = compressionType;
        this.inheritanceMap = inheritanceMap;

    }

    public void parse(ChildToParentMap childToParentMap) {

        haplotypeProbabilities = new HashMap<>(childToParentMap.children.length);
        missing = new HashSet<>();

        try {

            int tempInt = Integer.reverseBytes(compressedDataBlockContent.getInt());

            long tempLong = Integer.toUnsignedLong(tempInt);

            if (tempLong <= 0 || tempLong > Integer.MAX_VALUE) {

                throw new IllegalArgumentException("Unexpected uncompressed length (" + tempLong + ") for variant " + variantInformation.id + ", should be an integer between 0 (excluded) and " + Integer.MAX_VALUE + " (included).");

            }

            int uncompressedLength = (int) tempLong;

            byte[] compressedByteArray = new byte[blockLength - Integer.BYTES];

            compressedDataBlockContent.get(compressedByteArray, 0, blockLength - Integer.BYTES);

            byte[] uncompressedByteAray;

            if (compressionType == 2) {

                uncompressedByteAray = new byte[(int) uncompressedLength];

                long decompressedBytes = Zstd.decompress(uncompressedByteAray, compressedByteArray);

                if (decompressedBytes != uncompressedLength) {

                    throw new IllegalArgumentException(decompressedBytes + " bytes decompressed where " + uncompressedLength + " expected.");

                }

            } else if (compressionType == 1) {

                Inflater inflater = new Inflater(false);

                inflater.setInput(compressedByteArray);

                uncompressedByteAray = new byte[(int) uncompressedLength];

                int bytesUncompressed = inflater.inflate(uncompressedByteAray);

                if (bytesUncompressed == 0) {

                    throw new IllegalArgumentException("Missing input or dictionary.");

                } else if (bytesUncompressed != uncompressedLength) {

//                String debug = new String(uncompressedByteAray, 0, uncompressedByteAray.length, encoding);
                    throw new IllegalArgumentException("Unexpected number of bytes uncompressed " + bytesUncompressed + " (expected: " + uncompressedLength + ") for variant " + variantInformation.id + ".");

                }

            } else if (compressionType == 0) {

                uncompressedByteAray = compressedByteArray;

            } else {

                throw new IllegalArgumentException("Compression type " + compressionType + " not supported.");

            }

            ByteBuffer dataBlockContent = ByteBuffer.wrap(uncompressedByteAray);

            tempInt = Integer.reverseBytes(dataBlockContent.getInt());

            tempLong = Integer.toUnsignedLong(tempInt);

            if (tempLong <= 0 || tempLong > Integer.MAX_VALUE) {

                throw new IllegalArgumentException("Unexpected number of individuals (" + tempLong + ") for variant " + variantInformation.id + ", should be an integer between 0 (excluded) and " + Integer.MAX_VALUE + " (included).");

            }

            int nSamples = (int) tempLong;

            if (nSamples != sampleIds.length) {

                throw new IllegalArgumentException("Unexpected number of individuals (" + nSamples + ") found for variant " + variantInformation.id + ", " + sampleIds.length + " expected.");

            }

            short tempShort = Short.reverseBytes(dataBlockContent.getShort());
            nAlleles = Short.toUnsignedInt(tempShort);

            if (nAlleles != variantInformation.alleles.length) {

                throw new IllegalArgumentException("Unexpected number of alleles (" + nAlleles + ") found for variant " + variantInformation.id + ", " + variantInformation.alleles.length + " expected.");

            }

            int minPloidy = dataBlockContent.get();

            if (minPloidy < 0 || minPloidy > 63) {

                throw new IllegalArgumentException("Unexpected value for min ploidy (" + minPloidy + ") found for variant " + variantInformation.id + ", should be between " + 0 + " (inclusive) and " + 63 + " (inclusive).");

            }

            if (minPloidy <= 2 && nAlleles > 1) {

                int maxPloidy = dataBlockContent.get();

                if (maxPloidy < 0 || maxPloidy > 63) {

                    throw new IllegalArgumentException("Unexpected value for min ploidy (" + maxPloidy + ") found for variant " + variantInformation.id + ", should be between " + 0 + " (inclusive) and " + 63 + " (inclusive).");

                }

                int[] ploidyArray = new int[nSamples];

                for (int sampleI = 0; sampleI < nSamples; sampleI++) {

                    byte missingPloidy = dataBlockContent.get();

                    int missingValue = (missingPloidy >> 7) & 0b1;

                    if (missingValue != 0 && missingValue != 1) {

                        throw new IllegalArgumentException("Unexpected value for missingness (" + missingValue + ") found for variant " + variantInformation.id + " in sample " + sampleI + ", '0' or '1' expected.");

                    }

                    String sampleId = sampleIds[sampleI];

                    if (childToParentMap.sampleIds.contains(sampleId)) {

                        missing.add(sampleId);

                    }

                    int ploidyI = missingPloidy & 0b01111111;

                    if (ploidyI < minPloidy || ploidyI > maxPloidy) {

                        throw new IllegalArgumentException("Unexpected value for ploidy (" + ploidyI + ") found for variant " + variantInformation.id + " in sample " + sampleI + ", should be between " + minPloidy + " (inclusive) and " + maxPloidy + " (inclusive).");

                    }

                    ploidyArray[sampleI] = ploidyI;

                }

                int phased = dataBlockContent.get();

                if (phased != 0 && phased != 1) {

                    throw new IllegalArgumentException("Unexpected phased value (" + phased + ") found for variant " + variantInformation.id + ", '0' or '1' expected.");

                }

                if (phased == 0) {

                    throw new IllegalArgumentException("Genotypes not phased for variant " + variantInformation.id + ".");

                }

                int nBits = dataBlockContent.get();

                if (nBits < 1 || nBits > 32) {

                    throw new IllegalArgumentException("Unexpected value for the number of bits encoding the genotyping probabilities (" + nBits + ") found for variant " + variantInformation.id + ", should be between " + 1 + " (inclusive) and " + 32 + " (inclusive).");

                }

                int[] powers = new int[nBits];

                int value = 1;

                for (int i = 0; i < nBits; i++) {

                    powers[i] = value;

                    value *= 2;

                }

                double denominator = 2 * powers[nBits - 1] - 1;

                byte[] dataBlockBytes = new byte[dataBlockContent.remaining()];
                dataBlockContent.get(dataBlockBytes);

                BitSet bitSet = BitSet.valueOf(dataBlockBytes);

                int offSet = 0;

                for (int sampleI = 0; sampleI < nSamples; sampleI++) {

                    String sampleId = sampleIds[sampleI];

                    int z = ploidyArray[sampleI];

                    if (childToParentMap.sampleIds.contains(sampleId)) {

                        double[] probabilities = new double[(nAlleles - 1) * z];

                        int probabilityI = 0;

                        for (int contig = 1; contig <= z; contig++) {

                            for (int allele = 1; allele < nAlleles; allele++) {

                                double p = 0.0;

                                for (int i = 0; i < nBits; i++) {

                                    if (bitSet.get(offSet + probabilityI * nBits + nBits - i - 1)) {

                                        p += powers[i];

                                    }
                                }

                                p /= denominator;

                                probabilities[probabilityI] = p;

                                probabilityI++;

                            }
                        }

                        haplotypeProbabilities.put(sampleId, probabilities);

                    }

                    offSet += nBits * (nAlleles - 1) * z;

                }

                double[] nAltChildren = new double[variantInformation.alleles.length];
                double[] nContigChildren = new double[variantInformation.alleles.length];
                int nChildren = 0;

                for (String childId : childToParentMap.children) {

                    if (!missing.contains(childId)) {

                        nChildren++;

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            nAltChildren[alleleI] += getProbability(childId, alleleI);
                            nContigChildren[alleleI] += getPloidy(childId);

                        }
                    }
                }

                nChildrenGenotyped = nChildren;

                alleleFrequencyChildren = new double[variantInformation.alleles.length];

                for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                    alleleFrequencyChildren[alleleI] = nAltChildren[alleleI] / nContigChildren[alleleI];

                }
            }

        } catch (Exception e) {

            throw new RuntimeException(e);

        }
    }
    
    public boolean isMissing(String sampleId) {
        
        return missing.contains(sampleId);
        
    }

    public int getPloidy(
            String sampleId
    ) {

        double[] sampleProbabilities = haplotypeProbabilities.get(sampleId);
        return sampleProbabilities.length / (nAlleles - 1);

    }

    public double getProbability(
            String sampleId,
            int z,
            int alleleIndex
    ) {

        double[] sampleProbabilities = haplotypeProbabilities.get(sampleId);

        if (alleleIndex < nAlleles - 1) {

            return sampleProbabilities[z * (nAlleles - 1) + alleleIndex - 1];

        } else if (alleleIndex == nAlleles - 1) {

            double complement = 0.0;

            for (int otherAllele = 1; otherAllele < alleleIndex; otherAllele++) {

                complement += sampleProbabilities[z * (nAlleles - 1) + otherAllele - 1];

            }

            return Math.max(1.0 - complement, 0.0);

        } else {

            throw new IllegalArgumentException("Cannot retrieve allele index " + alleleIndex + " for sample " + sampleId + " and variant " + variantInformation.id + ", " + nAlleles + " alleles available.");

        }
    }

    public double getProbability(
            String sampleId,
            int alleleIndex
    ) {

        double p = 0.0;

        for (int z = 0; z < getPloidy(sampleId); z++) {

            p += getProbability(sampleId, z, alleleIndex);

        }

        return p;

    }

    public VariantInformation getVariantInformation() {

        return variantInformation;

    }

    public double[] getHaplotypes(
            String childId,
            String motherId,
            String fatherId,
            int testedAlleleIndex
    ) {

        if (!missing.contains(childId) && !(missing.contains(motherId) && missing.contains(fatherId))) {

            int ploidyChild = getPloidy(childId);

            char[] inheritance = inheritanceMap.get(ploidyChild);

            double motherTransmitted = 0.0;
            double fatherTransmitted = 0.0;

            for (int z = 0; z < ploidyChild; z++) {

                double probability = getProbability(fatherId, z, testedAlleleIndex);

                char parent = inheritance[z];

                switch (parent) {

                    case MOTHER:
                        motherTransmitted += probability;
                        break;

                    case FATHER:
                        fatherTransmitted += probability;
                        break;

                    default:
                        throw new IllegalArgumentException("Unsupported parent " + parent + ".");
                }
            }

            double mother = missing.contains(motherId) ? alleleFrequencyChildren[testedAlleleIndex] * getPloidy(motherId) : getProbability(motherId, testedAlleleIndex);

            double motherNonTransmitted = mother - motherTransmitted;

            double father = missing.contains(fatherId) ? alleleFrequencyChildren[testedAlleleIndex] * getPloidy(fatherId) : getProbability(fatherId, testedAlleleIndex);

            double fatherNonTransmitted = father - fatherTransmitted;

            return new double[]{motherNonTransmitted, motherTransmitted, fatherTransmitted, fatherNonTransmitted};

        }

        return null;

    }

    public double getAlleleFrequency(
            int testedAlleleIndex
    ) {
        return alleleFrequencyChildren[testedAlleleIndex];
    }

    public int getnChildrenGenotyped() {
        return nChildrenGenotyped;
    }

}
