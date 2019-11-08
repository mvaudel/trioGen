package no.uib.triogen.io.genotypes;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.io.genotypes.vcf.generic.VcfIterator;

/**
 * Enum of the implemented variant iterators.
 *
 * @author Marc Vaudel
 */
public enum GenotypesFileType {

    vcf(0, "VCF"),
    sangerVCF(1, "Sanger VCF");

    /**
     * The index of the file type.
     */
    public final int index;
    /**
     * The description of the file type.
     */
    public final String description;

    /**
     * Constructor.
     *
     * @param index the index of the file type
     * @param description the description of the file type
     */
    private GenotypesFileType(
            int index,
            String description
    ) {
        this.index = index;
        this.description = description;
    }

    /**
     * Returns the file type corresponding to the index.
     *
     * @param index the index
     *
     * @return the file type
     */
    public static GenotypesFileType getGenotypesFileType(
            int index
    ) {
        for (GenotypesFileType fileType : values()) {

            if (fileType.index == index) {

                return fileType;

            }
        }

        throw new IllegalArgumentException("No file type with index " + index + ".");

    }

    /**
     * Returns a variant iterator for the given file and file type.
     *
     * @param genotypesFile the file containing the genotypes
     * @param genotypesFileType the genotypes file type
     *
     * @return a variant iterator
     */
    public static VariantIterator getVariantIterator(
            File genotypesFile,
            GenotypesFileType genotypesFileType
    ) {

        switch (genotypesFileType) {
            case sangerVCF:
                return new CustomVcfIterator(genotypesFile);
            case vcf:
                return new VcfIterator(genotypesFile);
        }

        throw new IllegalArgumentException("No variant iterator implemented for " + genotypesFileType + " files.");

    }

    /**
     * Returns the different values as a string.
     * 
     * @return the different values as a string
     */
    public static String getCommandLineOptions() {

        return Arrays.stream(values())
                .map(
                        fileType -> String.join(", ", String.valueOf(fileType.index), fileType.description)
                )
                .collect(
                        Collectors.joining("; ")
                );

    }

}
