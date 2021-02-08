package no.uib.triogen.model.genome;

/**
 *
 *
 * @author Marc Vaudel
 */
public class VariantInformation {

    /**
     * The identifier of the variant.
     */
    public final String id;

    /**
     * The rsId of the variant.
     */
    public final String rsId;

    /**
     * The contig of the variant.
     */
    public final String contig;

    /**
     * The location of the variant.
     */
    public final int position;

    /**
     * The alleles.
     */
    public final String[] alleles;

    /**
     * The other alleles.
     */
    private String[] otherAlleles = null;

    /**
     * Constructor.
     *
     * @param id The identifier of the variant.
     * @param rsId The rsid of the variant.
     * @param contig The contig of the variant.
     * @param bp The location of the variant.
     * @param alleles The alleles.
     */
    public VariantInformation(
            String id,
            String rsId,
            String contig,
            int bp,
            String[] alleles
    ) {

        this.id = id;
        this.rsId = rsId;
        this.contig = contig;
        this.position = bp;
        this.alleles = alleles;

    }

    /**
     * Returns the other alleles as an array.
     * 
     * @param alleles The index of the allele of interest.
     * 
     * @return The other alleles as an array.
     */
    public static String[] getOtherAlleles(String[] alleles) {

        String[] otherAlleles = new String[alleles.length];

        for (int i = 0; i < alleles.length; i++) {

            StringBuilder otherAllele = new StringBuilder(2 * (alleles.length - 1) - 1);

            for (int j = 0; j < alleles.length; j++) {

                if (j != i) {

                    if (otherAllele.length() > 0) {

                        otherAllele.append(',');

                    }

                    otherAllele.append(alleles[j]);

                }
            }

            otherAlleles[i] = otherAllele.toString();

        }

        return otherAlleles;

    }

    /**
     * Returns a string containing the comma-separated list of other alleles for the given allele.
     * 
     * @param alleleIndex The index of the allele of interest.
     * 
     * @return A string containing the comma-separated list of other alleles for the given allele.
     */
    public String getOtherAllele(
            int alleleIndex
    ) {

        if (otherAlleles == null) {

            otherAlleles = getOtherAlleles(alleles);

        }

        return otherAlleles[alleleIndex];

    }
}
