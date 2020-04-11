package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Arrays;
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
    private float[] parentsP0sCache = null;
    /**
     * Cache for p0.
     */
    private double parentsP0Cache = Double.NaN;
    
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
            
            throw new IllegalArgumentException(altAlleles.size() + " alternative alleles found for variant " +  variantContext.getID() + ", only one supported.");
            
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
            
            throw new IllegalArgumentException(alleles.size() + " alleles found for variant " +  variantContext.getID() + " in sample " + sampleId + ", only one supported.");
            
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
        
        String attribute = variantContext.getAttributeAsString(DOSAGE_KEY, sampleId);
        
        if (attribute == null) {
            
            throw new IllegalArgumentException(
                    "Attribute for dosages (" + DOSAGE_KEY + ") not found."
            );
            
        }
        
        String[] split = attribute.split(",");
        
        float[] dosages = new float[3];
        
        for (int i = 0 ; i < 3 ; i++) {
            
            dosages[i] = Float.parseFloat(split[i]);
            
        }
        
        if (useCache) {
            
            dosagesCache.put(sampleId, dosages);
            
        }

        return dosages;
        
    }

    @Override
    public float[] getParentP0s(
            String[] childIds, 
            ChildToParentMap childToParentMap
    ) {

        if (parentsP0sCache == null) {

            float[] results = new float[2 * childIds.length];

            for (int i = 0; i < childIds.length; i++) {

                String childId = childIds[i];
                
                String motherId = childToParentMap.getMother(childId);
                results[i] = getDosages(motherId)[0];
                
                String fatherId = childToParentMap.getMother(childId);
                results[i + childIds.length] = getDosages(fatherId)[0];

            }
            
            parentsP0sCache = results;

        }
        
        return parentsP0sCache;

    }

    @Override
    public double getParentP0(String[] childIds, ChildToParentMap childToParentMap) {
        
        if (parentsP0Cache == Double.NaN) {

            double result = 0.0;

            for (float value : parentsP0sCache) {

                result += value;

            }
            
            parentsP0Cache = result;
            
        }
        
        return parentsP0Cache;
        
    }

}
