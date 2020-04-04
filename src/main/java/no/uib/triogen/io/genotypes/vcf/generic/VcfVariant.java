package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
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
     * The variant context.
     */
    private final VariantContext variantContext;
    /**
     * The alternative allele.
     */
    private Allele altAllele;
    
    /**
     * Constructor.
     * 
     * @param variantContext the variant context
     */
    public VcfVariant(VariantContext variantContext) {
        
        this.variantContext = variantContext;
        
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
    public short getGenotype(String sampleId) {
        
        Genotype genotype = variantContext.getGenotype(sampleId);
        
        List<Allele> alleles = genotype.getAlleles();
        
        if (alleles.size() != 2) {
            
            throw new IllegalArgumentException(alleles.size() + " alleles found for variant " +  variantContext.getID() + " in sample " + sampleId + ", only one supported.");
            
        }

        boolean allele11 = alleles.get(0).compareTo(altAllele) == 0;
        boolean allele21 = alleles.get(1).compareTo(altAllele) == 0;

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
    public String getVariantID() {
        
        return variantContext.getID();
        
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

}
