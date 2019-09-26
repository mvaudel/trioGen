package no.uib.triogen.io.vcf;

/**
 * This class contains and parses a line in the vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfLine {

    /**
     * The iterator used to parse the file.
     */
    private final VcfIterator vcfIterator;

    /**
     * The line as read from the file.
     */
    private final String line;

    /**
     * The indexes of the samples.
     */
    private int[] indexes;
    /**
     * The length of the variant description.
     */
    private int variantDescriptionLength = -1;

    /**
     * Constructor.
     *
     * @param vcfIterator the iterator used to parse the file
     * @param line the line as read from the file
     */
    public VcfLine(
            VcfIterator vcfIterator,
            String line
    ) {

        this.vcfIterator = vcfIterator;
        this.line = line;

    }

    /**
     * Parses the line.
     */
    public void parse() {

        indexes = new int[vcfIterator.getnSamples()];

        int nSeparators = 0;

        for (int index = 0; index < line.length(); index++) {

            if (line.charAt(index) == '\t') {

                nSeparators++;

                if (nSeparators >= vcfIterator.getnVariantColumns()) {

                    indexes[nSeparators - vcfIterator.getnVariantColumns()] = index;
                }
                if (nSeparators == 8) {

                    variantDescriptionLength = index;

                }
            }
        }
    }

    /**
     * Returns the description of the variant.
     *
     * @return the description of the variant
     */
    public String getVariantDescription() {

        return line.substring(0, variantDescriptionLength);

    }

    /**
     * Returns the genotype for a given sample. 0: 0|0 1: 1|0 2: 0|1 3: 1|1
     *
     * @param sampleId the id of the sample
     *
     * @return the genotype
     */
    public int getGenotype(
            String sampleId
    ) {

        int sampleIndex = vcfIterator.getSampleIndex(sampleId);
        int index1 = indexes[sampleIndex] + 1;
        int indexSeparator = index1 + 1;
        int index2 = indexSeparator + 1;
        char allele1 = line.charAt(index1);
        char separator = line.charAt(indexSeparator);
        char allele2 = line.charAt(index2);

        if (separator != '|' && separator != '/') {

            throw new IllegalArgumentException(
                    "Unexpected separator in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantDescription() + " in sample " + sampleId + "."
            );

        }
        if (allele1 != '0' && allele1 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allele1 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantDescription() + " in sample " + sampleId + "."
            );

        }
        if (allele2 != '0' && allele2 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allel2 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantDescription() + " in sample " + sampleId + "."
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
}
