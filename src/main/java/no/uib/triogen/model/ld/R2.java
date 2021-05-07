package no.uib.triogen.model.ld;

/**
 * This class contains the r2 value between the two alleles of two variants.
 *
 * @author Marc Vaudel
 */
public class R2 {

    public final int variantB;
    public final short alleleA;
    public final short alleleB;
    public final float r2Value;

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
}
