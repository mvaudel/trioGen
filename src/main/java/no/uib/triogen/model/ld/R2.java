package no.uib.triogen.model.ld;

/**
 * This class contains the r2 value between the two alleles of two variants.
 *
 * @author Marc Vaudel
 */
public class R2 {

    /**
     * The index of variant B in the LD matrix.
     */
    public final int variantB;
    /**
     * The id of variant B.
     */
    private String variantBId;
    /**
     * The rsid of variant B.
     */
    private String rsidBId;
    /**
     * The correlated allele for variant A.
     */
    public final short alleleA;
    /**
     * The correlated allele for variant B.
     */
    public final short alleleB;
    /**
     * The r2 value of the correlation.
     */
    public final float r2Value;

    /**
     * Constructor.
     * 
     * @param variantB The index of variant B in the LD matrix.
     * @param alleleA The correlated allele for variant A.
     * @param alleleB The correlated allele for variant B.
     * @param r2Value The r2 value of the correlation.
     */
    public R2(
            int variantB,
            short alleleA,
            short alleleB,
            float r2Value
    ) {

        this.variantB = variantB;
        this.alleleA = alleleA;
        this.alleleB = alleleB;
        this.r2Value = r2Value;

    }

    /**
     * Returns the id of variant B.
     * 
     * @return The id of variant B.
     */
    public String getVariantBId() {
        return variantBId;
    }

    /**
     * Sets the id of variant B.
     * 
     * @param variantBId The id of variant B.
     */
    public void setVariantBId(String variantBId) {
        this.variantBId = variantBId;
    }

    /**
     * Returns the rsid of variant B.
     * 
     * @return The esid of variant B.
     */
    public String getVariantBRsid() {
        return rsidBId;
    }

    /**
     * Sets the rsid of variant B.
     * 
     * @param rsidBId The esid of variant B.
     */
    public void setVariantBRsid(String rsidBId) {
        this.rsidBId = rsidBId;
    }
}
