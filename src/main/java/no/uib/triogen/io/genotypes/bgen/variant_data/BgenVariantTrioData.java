package no.uib.triogen.io.genotypes.bgen.variant_data;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.zip.Inflater;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import static no.uib.triogen.io.genotypes.InheritanceUtils.FATHER;
import static no.uib.triogen.io.genotypes.InheritanceUtils.MOTHER;
import no.uib.triogen.io.genotypes.bgen.BgenUtils;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.utils.CompressionUtils;

/**
 * Trio-based genotypes provider based on a bgen file.
 *
 * @author Marc Vaudel
 */
public class BgenVariantTrioData {

    /**
     * Array of the samples included.
     */
    private final String[] sampleIds;
    /**
     * Information on the variant.
     */
    private final VariantInformation variantInformation;
    /**
     * Byte buffer for the compressed content of the data block.
     */
    private final MappedByteBuffer compressedDataBlockContent;
    /**
     * The block length.
     */
    private final int blockLength;
    /**
     * Values for the haplotype probabilities for the children and parents.
     */
    private HashMap<String, double[]> haplotypeProbabilities = null;
    /**
     * Values for the allele probabilities for the children and parents summed across all contigs.
     */
    private HashMap<String, double[]> summedProbabilities = null;
    /**
     * The number of alleles.
     */
    private int nAlleles = -1;
    /**
     * The type of compression used.
     */
    private final int compressionType;
    /**
     * Estimate of the frequency of each allele.
     */
    private double[] alleleFrequency;
    /**
     * Indexes of the alleles ordered by decreasing frequency.
     */
    private int[] orderedAlleleIndex;
    /**
     * Boolean indicating whether alleles should be swapped for diploid
     * children.
     */
    private boolean swappedChildrenAllele = false;
    /**
     * The total number of children with genotyping information present.
     */
    private int nChildrenGenotyped;
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
     * @param sampleIds The sample ids of the bgen file.
     * @param variantInformation Information on the variant.
     * @param dataBlockContent The content of the data block.
     * @param blockLength The total block length.
     * @param compressionType The type of compression used.
     * @param inheritanceMap The allele inheritance map to use.
     * @param defaultMotherPloidy The default ploidy for mothers.
     * @param defaultFatherPloidy The default ploidy for fathers.
     */
    public BgenVariantTrioData(
            String[] sampleIds,
            VariantInformation variantInformation,
            MappedByteBuffer dataBlockContent,
            int blockLength,
            int compressionType,
            HashMap<Integer, char[]> inheritanceMap,
            int defaultMotherPloidy,
            int defaultFatherPloidy
    ) {

        this.sampleIds = sampleIds;
        this.variantInformation = variantInformation;
        this.compressedDataBlockContent = dataBlockContent;
        this.blockLength = blockLength;
        this.compressionType = compressionType;
        this.inheritanceMap = inheritanceMap;
        this.defaultMotherPloidy = defaultMotherPloidy;
        this.defaultFatherPloidy = defaultFatherPloidy;

    }

    /**
     * Parses the variant data.
     *
     * @param childToParentMap The child to parent map.
     * @param decompressor The decompressor to use.
     */
    public void parse(
            ChildToParentMap childToParentMap,
            ZstdDecompressor decompressor
    ) {

        haplotypeProbabilities = new HashMap<>(3 * childToParentMap.children.length);
        summedProbabilities = new HashMap<>(3 * childToParentMap.children.length);
        HashSet<String> missing = new HashSet<>();

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

                uncompressedByteAray = CompressionUtils.zstdDecompress(
                        decompressor, 
                        compressedByteArray, 
                        uncompressedLength
                );

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

                    if (missingValue == 1) {

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

                int[] powers = BgenUtils.getPowers(nBits);

                double denominator = 2 * powers[nBits - 1] - 1;

                byte[] dataBlockBytes = new byte[dataBlockContent.remaining()];
                dataBlockContent.get(dataBlockBytes);

                BitSet bitSet = BitSet.valueOf(dataBlockBytes);

                int offSet = 0;

                for (int sampleI = 0; sampleI < nSamples; sampleI++) {

                    String sampleId = sampleIds[sampleI];

                    int z = ploidyArray[sampleI];

                    if (!missing.contains(sampleId) || !childToParentMap.sampleIds.contains(sampleId)) {

                        double[] probabilities = new double[(nAlleles - 1) * z];
                        double[] sampleSummedProbabilities = new double[nAlleles];
                        double totalProbability = 0.0;

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

                                if (p < -0.1 || p > 1.1) {

                                    throw new IllegalArgumentException("Invalid probability (" + p + ") found for allele " + allele + " of contig " + contig + " in sample " + sampleId + ".");

                                }

                                probabilities[probabilityI] = p;
                                sampleSummedProbabilities[allele-1] += p;
                                totalProbability += p;

                                probabilityI++;

                            }
                        }
                        
                        sampleSummedProbabilities[nAlleles - 1] = z - totalProbability;

                        summedProbabilities.put(sampleId, sampleSummedProbabilities);
                        haplotypeProbabilities.put(sampleId, probabilities);

                    }

                    offSet += nBits * (nAlleles - 1) * z;

                }

                double[] nAltParents = new double[variantInformation.alleles.length];
                double[] nContigParents = new double[variantInformation.alleles.length];
                int nChildren = 0;

                for (String childId : childToParentMap.children) {

                    boolean childMissing = !haplotypeProbabilities.containsKey(childId);

                    if (!childMissing) {

                        nChildren++;

                    }

                    String motherId = childToParentMap.getMother(childId);
                    String fatherId = childToParentMap.getMother(childId);

                    boolean motherMissing = !haplotypeProbabilities.containsKey(motherId);
                    boolean fatherMissing = !haplotypeProbabilities.containsKey(fatherId);

                    if (!motherMissing) {

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            int ploidy = getPloidy(motherId);
                            nAltParents[alleleI] += getSummedProbability(motherId, alleleI);
                            nContigParents[alleleI] += ploidy;

                        }
                    }

                    if (!fatherMissing) {

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            int ploidy = getPloidy(fatherId);
                            nAltParents[alleleI] += getSummedProbability(fatherId, alleleI);
                            nContigParents[alleleI] += ploidy;

                        }
                    }

                    if (motherMissing && fatherMissing && !childMissing) {

                        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                            int ploidy = getPloidy(childId);
                            nAltParents[alleleI] += getSummedProbability(childId, alleleI);
                            nContigParents[alleleI] += ploidy;

                        }
                    }
                }

                nChildrenGenotyped = nChildren;

                alleleFrequency = new double[variantInformation.alleles.length];
                TreeMap<Double, ArrayList<Integer>> frequencyToAlleleIndexMap = new TreeMap<>();

                for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

                    double frequency = nAltParents[alleleI] / nContigParents[alleleI];
                    alleleFrequency[alleleI] = frequency;

                    ArrayList<Integer> allelesAtFrequency = frequencyToAlleleIndexMap.get(frequency);

                    if (allelesAtFrequency == null) {

                        allelesAtFrequency = new ArrayList<>(1);
                        frequencyToAlleleIndexMap.put(frequency, allelesAtFrequency);

                    }

                    allelesAtFrequency.add(alleleI);

                }

                orderedAlleleIndex = frequencyToAlleleIndexMap.descendingMap().entrySet().stream()
                        .flatMap(
                                entry -> entry.getValue().stream()
                        )
                        .mapToInt(
                                a -> a
                        )
                        .toArray();

            }
            
            IoUtils.closeBuffer(compressedDataBlockContent);

        } catch (Exception e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Returns a boolean indicating whether genotyping information is available
     * for the given sample.
     *
     * @param sampleId The id of the sample.
     *
     * @return A boolean indicating whether genotyping information is available
     * for the give sample is missing.
     */
    public boolean contains(String sampleId) {

        return haplotypeProbabilities.containsKey(sampleId);

    }

    /**
     * Returns the ploidy for the given sample.
     *
     * @param sampleId The id of the sample.
     *
     * @return The ploidy for the given sample.
     */
    public int getPloidy(
            String sampleId
    ) {

        double[] sampleProbabilities = haplotypeProbabilities.get(sampleId);
        return sampleProbabilities.length / (nAlleles - 1);

    }

    /**
     * Indicate that children alleles should be swapped.
     */
    public void swapChildrenAlleles() {

        swappedChildrenAllele = !swappedChildrenAllele;

    }

    /**
     * Returns a boolean indicating whether the given allele should be swapped
     * in children.
     *
     * @return A boolean indicating whether the given allele is swapped in
     * children.
     */
    public boolean isSwappedChildrenAlleles() {

        return swappedChildrenAllele;

    }

    /**
     * Returns the probability for the given allele in the given contig and
     * sample.
     *
     * @param sampleId The id of the sample.
     * @param z The index of the contig.
     * @param alleleIndex The index of the allele.
     *
     * @return The probability for the given allele in the given contig and
     * sample.
     */
    public double getProbability(
            String sampleId,
            int z,
            int alleleIndex
    ) {

        double[] sampleProbabilities = haplotypeProbabilities.get(sampleId);

        if (alleleIndex < nAlleles - 1) {

            return sampleProbabilities[z * (nAlleles - 1) + alleleIndex];

        } else if (alleleIndex == nAlleles - 1) {

            double complement = 0.0;

            for (int otherAllele = 0; otherAllele < alleleIndex; otherAllele++) {

                complement += sampleProbabilities[z * (nAlleles - 1) + otherAllele];

            }

            return Math.max(1.0 - complement, 0.0);

        } else {

            throw new IllegalArgumentException("Cannot retrieve allele index " + alleleIndex + " for sample " + sampleId + " and variant " + variantInformation.id + ", " + nAlleles + " alleles available.");

        }
    }

    /**
     * Returns the sum of probabilities for the given allele and sample over all
     * contigs.
     *
     * @param sampleId The id of the sample.
     * @param alleleIndex The index of the allele.
     *
     * @return The sum of probabilities for the given allele and sample over all
     * contigs.
     */
    public double getSummedProbability(
            String sampleId,
            int alleleIndex
    ) {

        double[] sampleProbabilities = summedProbabilities.get(sampleId);
        return sampleProbabilities[alleleIndex];

    }

    /**
     * Returns information on the variant.
     *
     * @return Information on the variant.
     */
    public VariantInformation getVariantInformation() {

        return variantInformation;

    }

    /**
     * Returns the haplotypes in an array: {motherNonTransmitted,
     * motherTransmitted, fatherTransmitted, fatherNonTransmitted}.Null if the
     * child or both parents are missing.If a parent is missing, haplotypes are
     * distributed according to the maf.
     *
     * @param childId The child id.
     * @param motherId The mother id.
     * @param fatherId The father id.
     * @param testedAlleleIndex The index of the allele tested in the variant
     * information.
     *
     * @return The haplotypes in an array.
     */
    public double[] getHaplotypes(
            String childId,
            String motherId,
            String fatherId,
            int testedAlleleIndex
    ) {

            int ploidyChild = getPloidy(childId);

            char[] inheritance = inheritanceMap.get(ploidyChild);

            double motherTransmitted = 0.0;
            double fatherTransmitted = 0.0;

            for (int z = 0; z < ploidyChild; z++) {

                double probability = getProbability(childId, z, testedAlleleIndex);

                char parent = inheritance[z];

                if (ploidyChild == 2 && swappedChildrenAllele) {

                    parent = InheritanceUtils.swap(parent);

                }

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

            double mother = !contains(motherId) ? alleleFrequency[testedAlleleIndex] * defaultMotherPloidy : getSummedProbability(motherId, testedAlleleIndex);

            double motherNonTransmitted = mother - motherTransmitted;

            double father = !contains(fatherId) ? alleleFrequency[testedAlleleIndex] * defaultFatherPloidy : getSummedProbability(fatherId, testedAlleleIndex);

            double fatherNonTransmitted = father - fatherTransmitted;

            return new double[]{motherNonTransmitted, motherTransmitted, fatherTransmitted, fatherNonTransmitted};

    }

    /**
     * Returns the estimate of the frequency of the given allele.
     *
     * @param testedAlleleIndex The index of the allele of interest.
     *
     * @return The estimate of the frequency of the given allele.
     */
    public double getAlleleFrequency(
            int testedAlleleIndex
    ) {
        return alleleFrequency[testedAlleleIndex];
    }

    /**
     * Returns the allele frequencies.
     * 
     * @return the allele frequencies
     */
    public double[] getAlleleFrequency() {
        return alleleFrequency;
    }

    /**
     * Returns the indexes of the alleles ordered by decreasing frequency.
     *
     * @return The indexes of the alleles ordered by decreasing frequency.
     */
    public int[] getOrderedAlleles() {

        return orderedAlleleIndex;

    }

    /**
     * Returns the number of genotyped children.
     *
     * @return the number of genotyped children
     */
    public int getnChildrenGenotyped() {
        return nChildrenGenotyped;
    }
}
