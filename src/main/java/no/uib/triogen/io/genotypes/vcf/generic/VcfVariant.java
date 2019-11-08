package no.uib.triogen.io.genotypes.vcf.generic;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.List;
import no.uib.triogen.io.genotypes.GenotypesProvider;

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
    public int getGenotype(String sampleId) {
        
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

}
