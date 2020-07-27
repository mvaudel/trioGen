package no.uib.triogen.io.genotypes.vcf.custom;

import java.util.HashMap;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * This class contains and parses a line in the vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfLine implements GenotypesProvider {

    /**
     * Vcf key for the genotyping flag.
     */
    public final static String TYPED = "TYPED";
    /**
     * The iterator used to parse the file.
     */
    private final CustomVcfIterator vcfIterator;

    /**
     * The line as read from the file.
     */
    private String line;
    /**
     * The variant id.
     */
    private String variantId;
    /**
     * The variant contig.
     */
    private String contig;
    /**
     * The variant location on the contig in bp.
     */
    private int bp;
    /**
     * The reference allele.
     */
    private String ref;
    /**
     * The alternative allele.
     */
    private String alt;
    /**
     * Boolean indicating whether the marker was genotyped.
     */
    private boolean typed;
    /**
     * Map of the indexes of samples in the alleles and dosages arrays.
     */
    private HashMap<String, Integer> indexMap;
    /**
     * Boolean array indicating whether the first allele is called alt.
     */
    private boolean[] alleles1;
    /**
     * Boolean array indicating whether the second allele is called alt.
     */
    private boolean[] alleles2;
    /**
     * Values for the dosages.
     */
    private float[][] dosages;
    /**
     * Cache for individual p0s.
     */
    private float[] parentsDosageP0sCache = null;
    /**
     * Cache for p0.
     */
    private double parentsDosageP0Cache = Double.NaN;
    /**
     * Cache for individual p0s.
     */
    private boolean[] parentsGenotypeP0sCache = null;
    /**
     * Cache for p0.
     */
    private int parentsGenotypeP0Cache = -1;

    /**
     * Constructor.
     *
     * @param vcfIterator The iterator used to parse the file.
     * @param line The line as read from the file.
     */
    public VcfLine(
            CustomVcfIterator vcfIterator,
            String line
    ) {

        this.vcfIterator = vcfIterator;
        this.line = line;

    }

    @Override
    public void parse(
            ChildToParentMap childToParentMap
    ) {

        indexMap = new HashMap<>(childToParentMap.sampleIds.size());
        alleles1 = new boolean[childToParentMap.sampleIds.size()];
        alleles2 = new boolean[childToParentMap.sampleIds.size()];
        dosages = new float[childToParentMap.sampleIds.size()][3];

        int nSeparators = 0;
        int previousIndex = -1;
        int vcfSampleIndex = 0;

        for (int index = 0; index < line.length(); index++) {

            if (line.charAt(index) == '\t') {

                nSeparators++;

                if (nSeparators == 1) {

                    contig = line.substring(previousIndex + 1, index);

                } else if (nSeparators == 2) {

                    bp = Integer.parseInt(line.substring(previousIndex + 1, index));

                } else if (nSeparators == 3) {

                    variantId = line.substring(previousIndex + 1, index);

                } else if (nSeparators == 4) {

                    ref = line.substring(previousIndex + 1, index);

                } else if (nSeparators == 5) {

                    alt = line.substring(previousIndex + 1, index);

                } else if (nSeparators == 8) {

                    String info = line.substring(previousIndex + 1, index);
                    typed = info.startsWith(TYPED);

                } else if (nSeparators > vcfIterator.getnVariantColumns()) {

                    String sampleId = vcfIterator.samples[vcfSampleIndex];

                    if (childToParentMap.sampleIds.contains(sampleId)) {

                        parseSample(
                                sampleId,
                                previousIndex,
                                index
                        );
                    }

                    vcfSampleIndex++;

                }

                previousIndex = index;

            }
        }

        String sampleId = vcfIterator.samples[vcfSampleIndex];

        if (childToParentMap.sampleIds.contains(sampleId)) {

            parseSample(
                    sampleId,
                    previousIndex,
                    line.length()
            );
        }

        line = null;

    }

    /**
     * Parses the values to keep in memory for the given sample.
     *
     * @param sampleId The id of the sample.
     * @param previousIndex The index of the separator before this sample.
     * @param index The index of the separator after this sample.
     */
    private void parseSample(
            String sampleId,
            int previousIndex,
            int index
    ) {

        int sampleIndex = indexMap.size();
        indexMap.put(sampleId, sampleIndex);

        // Parse genotypes 
        int indexAllele1 = previousIndex + 1;
        int indexSeparator = indexAllele1 + 1;
        int indexAllele2 = indexSeparator + 1;
        char allele1 = line.charAt(indexAllele1);
        char separator = line.charAt(indexSeparator);
        char allele2 = line.charAt(indexAllele2);

        if (separator != '|' && separator != '/' && separator != ':') {

            throw new IllegalArgumentException(
                    "Unexpected separator in genotype '" + allele1 + separator + allele2 + "' (variant '" + getVariantID() + "', line index '" + previousIndex + "')."
            );

        }
        if (allele1 != '0' && allele1 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allele1 in genotype '" + allele1 + separator + allele2 + "' (variant '" + getVariantID() + "', line index '" + previousIndex + "')."
            );

        }
        if (allele2 != '0' && allele2 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allel2 in genotype '" + allele1 + separator + allele2 + "' (variant '" + getVariantID() + "', line index '" + previousIndex + "')."
            );

        }

        alleles1[sampleIndex] = allele1 == '1';
        alleles2[sampleIndex] = allele2 == '1';

        // Parse dosages
        int indexBeginning = previousIndex + 1;
        int indexEnd = index;

        String sampleData = line.substring(indexBeginning, indexEnd);
        String[] split = sampleData.split(":");
        String dosagesString = split[split.length - 1];
        String[] dosagesStringSplit = dosagesString.split(",");

        if (dosagesStringSplit.length != 3) {

            throw new IllegalArgumentException(
                    dosagesStringSplit.length + " dosages found where 3 expected in '" + dosagesString + "' (variant '" + getVariantID() + "', line index '" + previousIndex + "')."
            );
        }

        float[] sampleDosages = new float[3];

        for (int i = 0; i < 3; i++) {

            try {

                sampleDosages[i] = Float.parseFloat(dosagesStringSplit[i]);

            } catch (Throwable t) {

                throw new IllegalArgumentException(
                        "Could not parse dosage '" + dosagesStringSplit[i] + "' as number in '" + dosagesString + "' (variant '" + getVariantID() + "', line index '" + previousIndex + "').",
                        t
                );
            }
        }

        dosages[sampleIndex] = sampleDosages;

    }

    @Override
    public String getVariantID() {

        return variantId;

    }

    @Override
    public String getContig() {

        return contig;

    }

    @Override
    public int getBp() {

        return bp;

    }

    @Override
    public String getRef() {

        return ref;

    }

    @Override
    public String getAlt() {

        return alt;

    }

    @Override
    public boolean genotyped() {

        return typed;

    }

    @Override
    public short getGenotype(
            String sampleId
    ) {

        int sampleIndex = indexMap.get(sampleId);

        boolean allele11 = alleles1[sampleIndex];
        boolean allele21 = alleles2[sampleIndex];

        if (!allele11 && !allele21) {

            return 0;

        } else if (allele11 && !allele21) {

            return 1;

        } else if (!allele11 && allele21) {

            return 2;

        } else {

            return 3;

        }
    }

    @Override
    public short getNAlt(
            String sampleId
    ) {

        int sampleIndex = indexMap.get(sampleId);

        boolean allele11 = alleles1[sampleIndex];
        boolean allele21 = alleles2[sampleIndex];

        if (!allele11 && !allele21) {

            return 0;

        } else if (allele11 && !allele21 || !allele11 && allele21) {

            return 1;

        } else {

            return 2;

        }
    }

    @Override
    public float[] getDosages(String sampleId) {

        int sampleIndex = indexMap.get(sampleId);

        return dosages[sampleIndex];

    }

    @Override
    public double getNAltDosages(
            String sampleId
    ) {

        float[] sampleDosages = getDosages(sampleId);

        return sampleDosages[1] + 2 * sampleDosages[2];

    }

    @Override
    public short[] getNAltH(
            String childId,
            String motherId,
            String fatherId
    ) {

        short genotypeChild = getGenotype(childId);
        short genotypeMother = getGenotype(motherId);
        short genotypeFather = getGenotype(fatherId);

        short nAltMother = (short) (genotypeMother >= 2 ? genotypeMother - 1 : genotypeMother);
        short nAltFather = (short) (genotypeFather >= 2 ? genotypeFather - 1 : genotypeFather);

        short h3 = (short) (genotypeChild == 0 || genotypeChild == 2 ? 0 : 1);
        short h1 = (short) (genotypeChild == 0 || genotypeChild == 1 ? 0 : 1);
        short h2 = (short) (nAltMother - h1);
        short h4 = (short) (nAltFather - h3);

        return new short[]{h1, h2, h3, h4};

    }

    @Override
    public void setParentP0s(
            String[] childIds,
            ChildToParentMap childToParentMap
    ) {

        float[] dosageHomRef = new float[2 * childIds.length];
        float[] dosageHomAlt = new float[2 * childIds.length];
        double sumDosageHomRef = 0.0;
        double sumDosageHomAlt = 0.0;
        boolean[] genotypeHomRef = new boolean[2 * childIds.length];
        boolean[] genotypeHomAlt = new boolean[2 * childIds.length];
        int sumGenotypeHomRef = 0;
        int sumGenotypeHomAlt = 0;

        for (int i = 0; i < childIds.length; i++) {

            String childId = childIds[i];

            String motherId = childToParentMap.getMother(childId);
            float[] sampleDosages = getDosages(motherId);
            dosageHomRef[i] = sampleDosages[0];
            dosageHomAlt[i] = sampleDosages[2];
            sumDosageHomRef += sampleDosages[0];
            sumDosageHomAlt += sampleDosages[2];

            int genotype = getGenotype(motherId);
            genotypeHomRef[i] = genotype == 0;
            genotypeHomAlt[i] = genotype == 3;

            if (genotype == 0) {

                sumGenotypeHomRef++;

            } else if (genotype == 3) {

                sumGenotypeHomAlt++;

            }

            String fatherId = childToParentMap.getFather(childId);
            sampleDosages = getDosages(fatherId);
            dosageHomRef[i + childIds.length] = sampleDosages[0];
            dosageHomAlt[i + childIds.length] = sampleDosages[2];
            sumDosageHomRef += sampleDosages[0];
            sumDosageHomAlt += sampleDosages[2];

            genotype = getGenotype(fatherId);
            genotypeHomRef[i + childIds.length] = genotype == 0;
            genotypeHomAlt[i + childIds.length] = genotype == 3;

            if (genotype == 0) {

                sumGenotypeHomRef++;

            } else if (genotype == 3) {

                sumGenotypeHomAlt++;

            }
        }

        parentsDosageP0sCache = sumDosageHomRef >= sumDosageHomAlt ? dosageHomRef : dosageHomAlt;
        parentsDosageP0Cache = sumDosageHomRef >= sumDosageHomAlt ? sumDosageHomRef : sumDosageHomAlt;

        parentsGenotypeP0sCache = sumGenotypeHomRef >= sumGenotypeHomAlt ? genotypeHomRef : genotypeHomAlt;
        parentsGenotypeP0Cache = sumGenotypeHomRef >= sumGenotypeHomAlt ? sumGenotypeHomRef : sumGenotypeHomAlt;

    }

    @Override
    public float[] getParentsDosageP0sCache() {

        return parentsDosageP0sCache;

    }

    @Override
    public double getParentsDosageP0Cache() {

        return parentsDosageP0Cache;

    }

    @Override
    public boolean[] getParentsGenotypeP0sCache() {

        return parentsGenotypeP0sCache;

    }

    @Override
    public int getParentsGenotypeP0Cache() {

        return parentsGenotypeP0Cache;

    }

    @Override
    public void emptyGenotypeDosageCaches() {

        indexMap = null;
        alleles1 = null;
        alleles2 = null;
        dosages = null;

    }
}
