package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.mendelian_error.MendelianErrorEstimator;

/**
 * The VcfVariant provides genotypes for a variant from a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfVariant implements GenotypesProvider {

    /**
     * Vcf key for the genotyping floag.
     */
    public final static String TYPED = "TYPED";
    /**
     * Vcf key for the dosages.
     */
    public final static String DOSAGE_KEY = "GP";

    /**
     * The variant context.
     */
    private final VariantContext variantContext;
    /**
     * The alternative allele.
     */
    private Allele altAllele;
    /**
     * Cache for the dosages.
     */
    private final HashMap<String, float[]> dosagesCache = new HashMap<>(0);
    /**
     * Cache for the hard calls.
     */
    private final HashMap<String, Short> hardCallsCache = new HashMap<>(0);
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
     * Constructor.
     *
     * @param variantContext The variant context.
     */
    public VcfVariant(
            VariantContext variantContext
    ) {

        this.variantContext = variantContext;

    }

    @Override
    public void parse(
            ChildToParentMap childToParentMap
    ) {

        // Get the alleles
        List<Allele> altAlleles = variantContext.getAlternateAlleles();

        if (altAlleles.size() != 1) {

            throw new IllegalArgumentException(altAlleles.size() + " alternative alleles found for variant " + variantContext.getID() + ", only one supported.");

        }

        altAllele = altAlleles.get(0);

        // Parse genotypes and dosages
        indexMap = new HashMap<>(childToParentMap.sampleIds.size());
        alleles1 = new boolean[childToParentMap.sampleIds.size()];
        alleles2 = new boolean[childToParentMap.sampleIds.size()];
        dosages = new float[childToParentMap.sampleIds.size()][3];

        for (String sampleId : childToParentMap.sampleIds) {

            int sampleIndex = indexMap.size();
            indexMap.put(sampleId, sampleIndex);

            Genotype genotype = variantContext.getGenotype(sampleId);

            // Genotypes 
            List<Allele> alleles = genotype.getAlleles();

            if (alleles.size() != 2) {

                throw new IllegalArgumentException(alleles.size() + " alleles found for variant " + variantContext.getID() + " in sample " + sampleId + ", only two supported.");

            }

            alleles1[sampleIndex] = alleles.get(0).compareTo(altAllele) == 0;
            alleles2[sampleIndex] = alleles.get(1).compareTo(altAllele) == 0;

            // Dosages
            Object dosagesObject = genotype.getAnyAttribute(DOSAGE_KEY);

            if (dosagesObject == null) {

                dosagesObject = genotype.getExtendedAttribute(DOSAGE_KEY);

                if (dosagesObject == null) {

                    dosagesObject = genotype.getLikelihoodsString();

                    if (dosagesObject == null) {

                        System.out.println("AD: " + genotype.getAD());
                        System.out.println("Filters: " + genotype.getFilters());
                        System.out.println("GenotypesString: " + genotype.getGenotypeString());
                        System.out.println("LikelihoodsString: " + genotype.getLikelihoodsString());
                        System.out.println("DP: " + genotype.getDP());
                        System.out.println("GQ: " + genotype.getGQ());
                        System.out.println("PL: " + genotype.getPL());
                        System.out.println("Ploidy: " + genotype.getPloidy());
                        System.out.println("Type: " + genotype.getType());

                        throw new IllegalArgumentException("No likelihood found for variant " + variantContext.getID() + " in sample " + sampleId + ".");

                    }
                }
            }

            String dosagesString = dosagesObject.toString();
            String[] dosagesSplit = dosagesString.split(",");

            if (dosagesSplit.length != 3) {

                throw new IllegalArgumentException(
                        dosagesSplit.length + " dosages found where 3 expected in '" + dosagesString + "' (variant '" + getVariantID() + "')."
                );
            }

            float[] sampleDosages = new float[3];

            for (int i = 0; i < 3; i++) {

                try {

                    sampleDosages[i] = Float.parseFloat(dosagesSplit[i]);

                } catch (Throwable t) {

                    throw new IllegalArgumentException(
                            "Could not parse dosage '" + dosagesSplit[i] + "' as number in '" + dosagesString + "' (variant '" + getVariantID() + "').",
                            t
                    );
                }
            }

            dosages[sampleIndex] = sampleDosages;

        }
    }

    @Override
    public String getVariantID() {

        return variantContext.getID();

    }

    @Override
    public String getContig() {

        return variantContext.getContig();

    }

    @Override
    public int getBp() {

        return variantContext.getStart();

    }

    @Override
    public String getRef() {

        return variantContext.getReference().getBaseString();

    }

    @Override
    public String getAlt() {

        return variantContext.getAlternateAlleles().stream()
                .map(
                        allele -> allele.getBaseString()
                )
                .collect(
                        Collectors.joining(",")
                );
    }

    @Override
    public boolean genotyped() {

        return variantContext.hasAttribute(TYPED);

    }

    @Override
    public short getGenotype(String sampleId) {

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

        short genotypeKid = getGenotype(childId);
        short genotypeMother = getGenotype(motherId);
        short genotypeFather = getGenotype(fatherId);

        short nAltMother = (short) (genotypeMother >= 2 ? genotypeMother - 1 : genotypeMother);
        short nAltFather = (short) (genotypeFather >= 2 ? genotypeFather - 1 : genotypeFather);

        short h3 = (short) (genotypeKid == 0 || genotypeKid == 2 ? 0 : 1);
        short h1 = (short) (genotypeKid == 0 || genotypeKid == 1 ? 0 : 1);
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
            dosageHomRef[i + childIds.length] = dosages[0];
            dosageHomAlt[i + childIds.length] = dosages[2];
            sumDosageHomRef += dosages[0];
            sumDosageHomAlt += dosages[2];

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

    @Override
    public double checkMendelianErrors(ChildToParentMap childToParentMap) {
        
        return MendelianErrorEstimator.estimateMendelianErrorPrevalence(this, childToParentMap);
        
    }

}
