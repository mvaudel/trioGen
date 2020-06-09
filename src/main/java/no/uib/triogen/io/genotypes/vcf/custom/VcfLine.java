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
     * Vcf key for the genotyping floag.
     */
    public final static String TYPED = "TYPED";
    /**
     * The iterator used to parse the file.
     */
    private final CustomVcfIterator vcfIterator;

    /**
     * The line as read from the file.
     */
    private final String line;

    /**
     * The indexes of the samples.
     */
    private int[] indexes;
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
     * Boolean indicating whether the genotypes provider should cache genotype
     * values.
     */
    private final boolean useCache;
    /**
     * Cache for the dosages.
     */
    private final HashMap<String, float[]> dosagesCache = new HashMap<>(0);
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
     * @param useCache Boolean indicating whether the genotypes provider should
     * cache genotype values.
     */
    public VcfLine(
            CustomVcfIterator vcfIterator,
            String line,
            boolean useCache
    ) {

        this.vcfIterator = vcfIterator;
        this.line = line;
        this.useCache = useCache;

    }

    @Override
    public void parse() {

        indexes = new int[vcfIterator.getnSamples()];

        int nSeparators = 0;
        int index1 = -1;

        for (int index = 0; index < line.length(); index++) {

            if (line.charAt(index) == '\t') {

                nSeparators++;

                if (nSeparators >= vcfIterator.getnVariantColumns()) {

                    indexes[nSeparators - vcfIterator.getnVariantColumns()] = index;
                }
                if (nSeparators == 1) {

                    contig = line.substring(index1 + 1, index);
                    index1 = index;

                } else if (nSeparators == 2) {

                    bp = Integer.parseInt(
                            line.substring(index1 + 1, index)
                    );
                    index1 = index;

                } else if (nSeparators == 3) {

                    variantId = line.substring(index1 + 1, index);
                    index1 = index;

                } else if (nSeparators == 4) {

                    ref = line.substring(index1 + 1, index);
                    index1 = index;

                } else if (nSeparators == 5) {

                    alt = line.substring(index1 + 1, index);
                    index1 = index;

                } else if (nSeparators == 7) {

                    index1 = index;

                } else if (nSeparators == 8) {

                    String info = line.substring(index1 + 1, index);
                    typed = info.startsWith(TYPED);
                    index1 = index;

                }
            }
        }
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

        int sampleIndex = vcfIterator.getSampleIndex(sampleId);
        int index1 = indexes[sampleIndex] + 1;
        int indexSeparator = index1 + 1;
        int index2 = indexSeparator + 1;
        char allele1 = line.charAt(index1);
        char separator = line.charAt(indexSeparator);
        char allele2 = line.charAt(index2);

        if (separator != '|' && separator != '/' && separator != ':') {

            throw new IllegalArgumentException(
                    "Unexpected separator in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID() + " in sample " + sampleId + "."
            );

        }
        if (allele1 != '0' && allele1 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allele1 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID() + " in sample " + sampleId + "."
            );

        }
        if (allele2 != '0' && allele2 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allel2 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID() + " in sample " + sampleId + "."
            );

        }

        boolean allele11 = allele1 == '1';
        boolean allele21 = allele2 == '1';

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
    public short[] getH(
            ChildToParentMap childToParentMap,
            String childId
    ) {

        String motherId = childToParentMap.getMother(childId);
        String fatherId = childToParentMap.getFather(childId);
        
        return getH(
                childId, 
                motherId, 
                fatherId
        );

    }

    @Override
    public short[] getH(
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
    public float[] getDosages(String sampleId) {

        if (useCache) {

            float[] dosages = dosagesCache.get(sampleId);

            if (dosages != null) {

                return dosages;

            }
        }

        int sampleIndex = vcfIterator.getSampleIndex(sampleId);
        int index1 = indexes[sampleIndex] + 1;
        int index2 = sampleIndex == indexes.length - 1 ? line.length() : indexes[sampleIndex + 1];

        String sampleData = line.substring(index1, index2);
        String[] split = sampleData.split(":");
        String[] dosagesString = split[split.length - 1].split(",");

        if (dosagesString.length != 3) {

            throw new IllegalArgumentException(
                    dosagesString.length + " dosages found where 3 expected."
            );
        }

        float[] dosages = new float[3];

        for (int i = 0; i < 3; i++) {

            try {

                dosages[i] = Float.parseFloat(dosagesString[i]);

            } catch (Throwable t) {

                throw new IllegalArgumentException(
                        "Could not parse dosage as number (Sample: '" + sampleId + "' value: '" + dosagesString[i] + "')",
                        t
                );
            }
        }

        if (useCache) {

            dosagesCache.put(sampleId, dosages);

        }

        return dosages;

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
            float[] dosages = getDosages(motherId);
            dosageHomRef[i] = dosages[0];
            dosageHomAlt[i] = dosages[2];
            sumDosageHomRef += dosages[0];
            sumDosageHomAlt += dosages[2];
            
            int genotype = getGenotype(motherId);
            genotypeHomRef[i] = genotype == 0;
            genotypeHomAlt[i] = genotype == 3;

            if (genotype == 0) {
                
                sumGenotypeHomRef++;
            
            } else if (genotype == 3) {
                
                sumGenotypeHomAlt++;
            
            }

            String fatherId = childToParentMap.getFather(childId);
            dosages = getDosages(fatherId);
            dosageHomRef[i] = dosages[0];
            dosageHomAlt[i] = dosages[2];
            sumDosageHomRef += dosages[0];
            sumDosageHomAlt += dosages[2];
            
            genotype = getGenotype(fatherId);
            genotypeHomRef[i] = genotype == 0;
            genotypeHomAlt[i] = genotype == 3;

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
}
