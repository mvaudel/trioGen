package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.GenotypesFileType;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * This script test that two parsers give the same results.
 *
 * @author Marc Vaudel
 */
public class TestVcfParsers {

    public final static int N_VARIANTS = 1000000;
    public final static double DOSAGE_TOLERANCE = 1e-6;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        if (args.length != 2) {

            throw new IllegalArgumentException("Two arguments expected: vcf_file trio_file");

        }

        File vcfFile = new File(args[0]);

        if (!vcfFile.exists()) {

            throw new IllegalArgumentException("VCF file " + vcfFile.getAbsolutePath() + " not found.");

        }

        File trioFile = new File(args[1]);

        if (!trioFile.exists()) {

            throw new IllegalArgumentException("Trio file " + trioFile.getAbsolutePath() + " not found.");

        }

        try {

            ChildToParentMap childToParentMap = ChildToParentMap.fromFile(trioFile);

            int nVariants = 0;

            VariantIterator genericIterator = GenotypesFileType.getVariantIterator(
                    vcfFile,
                    GenotypesFileType.vcf,
                    null,
                    0
            );
            VariantIterator customIterator = GenotypesFileType.getVariantIterator(
                    vcfFile,
                    GenotypesFileType.sangerVCF,
                    null,
                    0
            );

            try {

                int lastProgress = 0;

                while (nVariants++ < N_VARIANTS) {

                    int currentProgress = 100 * nVariants / N_VARIANTS;

                    if (currentProgress > lastProgress) {

                        lastProgress = currentProgress;

                        System.out.println(Instant.now() + "    " + nVariants + " processed (" + currentProgress + "%)");

                    }

                    GenotypesProvider genericProvider = genericIterator.next();
                    genericProvider.parse(childToParentMap);

                    GenotypesProvider customProvider = customIterator.next();
                    customProvider.parse(childToParentMap);

                    if (!genericProvider.getVariantID().equals(customProvider.getVariantID())) {

                        throw new IllegalArgumentException(
                                "Different variants at line " + nVariants + ":\n"
                                + "Generic: " + genericProvider.getVariantID() + ":\n"
                                + "Custom: " + customProvider.getVariantID()
                        );

                    }

                    if (!genericProvider.getRef().equals(customProvider.getRef())) {

                        throw new IllegalArgumentException(
                                "Ref allele mismatch at line " + nVariants + ":\n"
                                + "Generic: " + genericProvider.getRef() + ":\n"
                                + "Custom: " + customProvider.getRef()
                        );

                    }

                    if (!genericProvider.getAlt().equals(customProvider.getAlt())) {

                        throw new IllegalArgumentException(
                                "Alt allele mismatch at line " + nVariants + ":\n"
                                + "Generic: " + genericProvider.getAlt() + ":\n"
                                + "Custom: " + customProvider.getAlt()
                        );

                    }

                    if (genericProvider.genotyped() != customProvider.genotyped()) {

                        throw new IllegalArgumentException(
                                "Typing mismatch at line " + nVariants + ":\n"
                                + "Generic: " + genericProvider.genotyped() + ":\n"
                                + "Custom: " + customProvider.genotyped()
                        );

                    }

                    for (String id : childToParentMap.sampleIds) {

                        if (genericProvider.getNAlt(id) != customProvider.getNAlt(id)) {

                            throw new IllegalArgumentException(
                                    "nAlt mismatch for sample " + id + " at line " + nVariants + ":\n"
                                    + "Generic: " + genericProvider.getNAlt(id) + ":\n"
                                    + "Custom: " + customProvider.getNAlt(id)
                            );

                        }

                        if (genericProvider.getGenotype(id) != customProvider.getGenotype(id)) {

                            throw new IllegalArgumentException(
                                    "Genotyping mismatch for sample " + id + " at line " + nVariants + ":\n"
                                    + "Generic: " + genericProvider.getGenotype(id) + ":\n"
                                    + "Custom: " + customProvider.getGenotype(id)
                            );

                        }

                        for (int i = 0; i < 3; i++) {

                            double genericDosage = genericProvider.getDosages(id)[i];
                            double customDosage = customProvider.getDosages(id)[i];

                            if (!Double.isFinite(genericDosage) || !Double.isFinite(customDosage) || Math.abs(customDosage - genericDosage) > DOSAGE_TOLERANCE) {

                                throw new IllegalArgumentException(
                                        "Dosage mismatch for sample " + id + " at line " + nVariants + ":\n"
                                        + "Generic: " + genericProvider.getDosages(id)[i] + ":\n"
                                        + "Custom: " + customProvider.getDosages(id)[i]
                                );

                            }
                        }
                    }
                }

            } finally {

                genericIterator.close();
                customIterator.close();

            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }
}
