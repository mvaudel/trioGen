package no.uib.triogen.io.genotypes.vcf.reader;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.mendelian_error.MendelianErrorEstimator;

/**
 * The VcfVariant provides genotypes for a variant from a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfVariant {

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
    
    /**
     * Returns information on this variant.
     * 
     * @return information on this variant
     */
    public VariantInformation getVariantInformation() {
        
        String contig = variantContext.getContig();
        
        int bp = variantContext.getStart();
        
        String[] alleles = variantContext.getAlleles().stream()
                .map(
                allele -> allele.getBaseString()
                )
                .toArray(String[]::new);
        
        String id = String.join("_", contig, Integer.toString(bp), String.join("_", alleles));
        
        String rsId = variantContext.getID().startsWith("rs") && !variantContext.getID().contains("_") ? variantContext.getID() : "";
        
        return new VariantInformation(id, rsId, contig, bp, alleles);
        
    }
    
    /**
     * Returns the phased genotypes for the given samples.
     * 
     * @param sampleIds the sample ids
     * 
     * @return the phased genotypes
     */
    public ArrayList<String[]> getGenotypes(String[] sampleIds) {
        
        ArrayList<String[]> result = new ArrayList<>(sampleIds.length);
        
        for (String sampleId : sampleIds) {

            Genotype genotype = variantContext.getGenotype(sampleId);

            List<Allele> alleles = genotype.getAlleles();
            
            String[] sampleAlleles = new String[alleles.size()];
            
            for (int i = 0 ; i < alleles.size() ; i++) {
                
                sampleAlleles[i] = alleles.get(i).getBaseString();
                
            }
            
            result.add(sampleAlleles);
            
        }
        
        return result;
        
    }

    /**
     * Parses genotypes and dosages for the given samples.
     * 
     * @param childToParentMap the child to parents map
     */
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

            alleles1[sampleIndex] = alleles.get(1).compareTo(altAllele) == 0;
            alleles2[sampleIndex] = alleles.get(0).compareTo(altAllele) == 0;

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
                        dosagesSplit.length + " dosages found where 3 expected in '" + dosagesString + "' (variant '" + variantContext.getID() + "')."
                );
            }

            float[] sampleDosages = new float[3];

            for (int i = 0; i < 3; i++) {

                try {

                    sampleDosages[i] = Float.parseFloat(dosagesSplit[i]);

                } catch (Throwable t) {

                    throw new IllegalArgumentException(
                            "Could not parse dosage '" + dosagesSplit[i] + "' as number in '" + dosagesString + "' (variant '" + variantContext.getID() + "').",
                            t
                    );
                }
            }

            dosages[sampleIndex] = sampleDosages;

        }
    }

    /**
     * Returns a boolean indicating whether the variant is genotyped.
     * 
     * @return a boolean indicating whether the variant is genotyped
     */
    public boolean genotyped() {

        return variantContext.hasAttribute(TYPED);

    }

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

    public float[] getGenotypingProbabilities(String sampleId) {

        int sampleIndex = indexMap.get(sampleId);

        return dosages[sampleIndex];

    }

    public double getNAltGenotypingProbabilities(
            String sampleId
    ) {

        float[] sampleDosages = getGenotypingProbabilities(sampleId);

        return sampleDosages[1] + 2 * sampleDosages[2];

    }

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
            float[] dosages = getGenotypingProbabilities(motherId);
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
            dosages = getGenotypingProbabilities(fatherId);
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

    public float[] getParentsDosageP0sCache() {

        return parentsDosageP0sCache;

    }

    public double getParentsDosageP0Cache() {

        return parentsDosageP0Cache;

    }

    public boolean[] getParentsGenotypeP0sCache() {

        return parentsGenotypeP0sCache;

    }

    public int getParentsGenotypeP0Cache() {

        return parentsGenotypeP0Cache;

    }

    public void emptyGenotypeDosageCaches() {

        indexMap = null;
        alleles1 = null;
        alleles2 = null;
        dosages = null;

    }

    public String getGenotypingProbabilitiesAsString(
            String sampleId,
            String separator
    ) {

        float[] dosages = getGenotypingProbabilities(sampleId);

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < dosages.length; i++) {

            if (i > 0) {

                sb.append(separator);

            }

            sb.append(dosages[i]);

        }

        return sb.toString();

    }
}
