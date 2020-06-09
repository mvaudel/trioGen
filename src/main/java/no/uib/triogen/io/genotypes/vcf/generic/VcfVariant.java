package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;

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
     * Boolean indicating whether the genotypes provider should cache genotype
     * values.
     */
    private final boolean useCache;
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
     * Constructor.
     *
     * @param variantContext The variant context.
     * @param useCache Boolean indicating whether the genotypes provider should
     * cache genotype values.
     */
    public VcfVariant(
            VariantContext variantContext,
            boolean useCache
    ) {

        this.variantContext = variantContext;
        this.useCache = useCache;

    }

    @Override
    public void parse() {

        List<Allele> altAlleles = variantContext.getAlternateAlleles();

        if (altAlleles.size() != 1) {

            throw new IllegalArgumentException(altAlleles.size() + " alternative alleles found for variant " + variantContext.getID() + ", only one supported.");

        }

        altAllele = altAlleles.get(0);

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

        if (useCache) {

            Short genotype = hardCallsCache.get(sampleId);

            if (genotype != null) {

                return genotype;

            }
        }

        Genotype genotype = variantContext.getGenotype(sampleId);

        List<Allele> alleles = genotype.getAlleles();

        if (alleles.size() != 2) {

            throw new IllegalArgumentException(alleles.size() + " alleles found for variant " + variantContext.getID() + " in sample " + sampleId + ", only two supported.");

        }

        boolean allele11 = alleles.get(0).compareTo(altAllele) == 0;
        boolean allele21 = alleles.get(1).compareTo(altAllele) == 0;

        if (!allele11 && !allele21) {

            short result = 0;

            if (useCache) {

                hardCallsCache.put(sampleId, result);

            }

            return result;

        } else if (allele11 && !allele21) {

            short result = 1;

            if (useCache) {

                hardCallsCache.put(sampleId, result);

            }

            return result;

        } else if (!allele11 && allele21) {

            short result = 2;

            if (useCache) {

                hardCallsCache.put(sampleId, result);

            }

            return result;

        } else {

            short result = 3;

            if (useCache) {

                hardCallsCache.put(sampleId, result);

            }

            return result;

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
    public float[] getDosages(String sampleId) {

        if (useCache) {

            float[] dosages = dosagesCache.get(sampleId);

            if (dosages != null) {

                return dosages;

            }
        }

        Genotype genotype = variantContext.getGenotype(sampleId);

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
        
        String[] dosagesSplit = dosagesObject.toString().split(",");

        if (dosagesSplit.length != 3) {

            throw new IllegalArgumentException(
                    dosagesSplit.length + " dosages found where 3 expected."
            );
        }
        
        float[] dosages = new float[3];

        for (int i = 0; i < 3; i++) {

            try {

                dosages[i] = Float.parseFloat(dosagesSplit[i]);

            } catch (Throwable t) {

                throw new IllegalArgumentException(
                        "Could not parse dosage as number (Sample: '" + sampleId + "' value: '" + dosagesSplit[i] + "')",
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
