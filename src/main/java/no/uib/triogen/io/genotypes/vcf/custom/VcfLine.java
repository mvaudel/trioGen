package no.uib.triogen.io.genotypes.vcf.custom;

import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * This class contains and parses a line in the vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfLine implements GenotypesProvider {

    /**
     * The iterator used to parse the file.
     */
    private final CustomVcfIterator vcfIterator;

    /**
     * The line as read from the file.
     */
    private final String line;

    /**
     * The indexes of the samples.
     */
    private int[] indexes;
    /**
     * The variant id.
     */
    private String variantId;

    /**
     * Constructor.
     *
     * @param vcfIterator the iterator used to parse the file
     * @param line the line as read from the file
     */
    public VcfLine(
            CustomVcfIterator vcfIterator,
            String line
    ) {

        this.vcfIterator = vcfIterator;
        this.line = line;

    }

    @Override
    public void parse() {

        indexes = new int[vcfIterator.getnSamples()];

        int nSeparators = 0;
        int index1 = 0;

        for (int index = 0; index < line.length(); index++) {

            if (line.charAt(index) == '\t') {

                nSeparators++;

                if (nSeparators >= vcfIterator.getnVariantColumns()) {

                    indexes[nSeparators - vcfIterator.getnVariantColumns()] = index;
                }
                if (nSeparators == 2) {

                    index1 = index;

                } else if (nSeparators == 3) {
                    
                    variantId = line.substring(index1 + 1, index);
                
                }
            }
        }
    }

    @Override
    public String getVariantID() {

        return variantId;

    }

    @Override
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
                    "Unexpected separator in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID()+ " in sample " + sampleId + "."
            );

        }
        if (allele1 != '0' && allele1 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allele1 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID() + " in sample " + sampleId + "."
            );

        }
        if (allele2 != '0' && allele2 != '1') {

            throw new IllegalArgumentException(
                    "Unexpected allel2 in genotype " + line.substring(index1, index2 + 1) + " for variant " + getVariantID() + " in sample " + sampleId + "."
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
    
    @Override
    public double[] getH(
            ChildToParentMap childToParentMap,
            String childId
    ) {

        String motherId = childToParentMap.getMother(childId);
        String fatherId = childToParentMap.getFather(childId);

        int genotypeChild = getGenotype(childId);
        int genotypeMother = getGenotype(motherId);
        int genotypeFather = getGenotype(fatherId);

        int nAltMother = genotypeMother >= 2 ? genotypeMother - 1 : genotypeMother;
        int nAltFather = genotypeFather >= 2 ? genotypeFather - 1 : genotypeFather;

        double h3 = genotypeChild == 0 || genotypeChild == 2 ? 0.0 : 1.0;
        double h1 = genotypeChild == 0 || genotypeChild == 1 ? 0.0 : 1.0;
        double h2 = nAltMother - h1;
        double h4 = nAltFather - h3;

        return new double[]{h1, h2, h3, h4};

    }
}
