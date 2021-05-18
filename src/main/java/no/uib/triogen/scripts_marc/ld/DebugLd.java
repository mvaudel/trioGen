package no.uib.triogen.scripts_marc.ld;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.processing.ld.P0Cache;

/**
 * Debugs the LD calculation.
 *
 * @author Marc Vaudel
 */
public class DebugLd {

    /**
     * The allele frequency threshold to use.
     */
    private static final double alleleFrequencyThreshold = 0.001;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            String chromosome = "20";
            String rsid = "rs6040450";
            String variantId = "20:11209782_C_T";
            int position = 11209782;

            File output = new File("/mnt/work/marc/moba/trioGen/tmp");

            File trioFile = new File("/mnt/work/marc/moba/run/triogen/pheno/trio");

            File bgenFile = new File("/mnt/archive2/marc/moba/phased_bgen/20.phased.bgen");
            HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(chromosome);
            int defaultMotherPloidy = InheritanceUtils.getDefaultMotherPloidy(chromosome);
            int defaultFatherPloidy = InheritanceUtils.getDefaultFatherPloidy(chromosome);

            int maxDistance = 500000;

            P0Cache p0Cache = new P0Cache(1);

            System.out.println("Parsing " + bgenFile.getAbsolutePath() + ".");

            long start = Instant.now().getEpochSecond();

            BgenIndex bgenIndex = BgenIndex.getBgenIndex(bgenFile);

            BgenFileReader bgenFileReader = new BgenFileReader(
                    bgenFile,
                    bgenIndex,
                    inheritanceMap,
                    defaultMotherPloidy,
                    defaultFatherPloidy
            );

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            System.out.println("Parsing " + bgenFile + " done (" + duration + " seconds)");

            System.out.println("Parsing " + trioFile.getAbsolutePath() + ".");

            start = Instant.now().getEpochSecond();

            ChildToParentMap childToParentMap = ChildToParentMap.fromFile(trioFile);

            end = Instant.now().getEpochSecond();
            duration = end - start;

            System.out.println("Parsing " + trioFile + " done (" + duration + " seconds)");

            System.out.println("Loading results for " + variantId + ".");

            start = Instant.now().getEpochSecond();

            ZstdDecompressor decompressor = new ZstdDecompressor();

            int indexA = -1;
            VariantInformation variantInformationA = null;

            for (int i = 0; i < bgenIndex.variantIdArray.length; i++) {

                variantInformationA = bgenIndex.variantInformationArray[i];

                if (variantInformationA.rsId.equals(rsid) || variantInformationA.id.equals(variantId)) {

                    indexA = i;

                    break;

                }
            }

            if (indexA == -1) {

                System.out.println("Variant not found. Variants within 10kb");

                for (int i = 0; i < bgenIndex.variantIdArray.length; i++) {

                    variantInformationA = bgenIndex.variantInformationArray[i];

                    if (variantInformationA.position >= position - 10000 && variantInformationA.position <= position + 10000) {

                        System.out.println(variantInformationA.id);

                    }
                }

                throw new IllegalArgumentException("Variant " + rsid + " not found.");

            }

            float[][] pHomA = p0Cache.getPHomozygous(variantInformationA.id);

            if (pHomA == null) {

                BgenVariantData variantData = bgenFileReader.getVariantData(indexA);
                variantData.parse(
                        childToParentMap,
                        decompressor
                );

                if (!hasAlleles(variantData)) {

                    throw new IllegalArgumentException("No allele found for " + rsid + ".");

                }

                p0Cache.register(variantData, childToParentMap);

            }

            ArrayList<R2> r2s = new ArrayList<>(2);

            pHomA = p0Cache.getPHomozygous(variantInformationA.id);

            try (SimpleFileWriter writer = new SimpleFileWriter(output, false)) {

                writer.writeLine(
                        "variant_A",
                        "rsid_A",
                        "allele_A",
                        "variant_B",
                        "rsid_A",
                        "allele_B",
                        "d",
                        "r2"
                );

                VariantIterator iteratorB = new VariantIterator(bgenIndex, variantInformationA.position - maxDistance, variantInformationA.position + maxDistance);

                Integer indexB;
                while ((indexB = iteratorB.next()) != null) {

                    VariantInformation variantInformationB = bgenIndex.variantInformationArray[indexB];

                    float[][] pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                    if (pHomB == null) {

                        BgenVariantData variantData = bgenFileReader.getVariantData(indexB);
                        variantData.parse(
                                childToParentMap,
                                decompressor
                        );

                        if (!hasAlleles(variantData)) {

                            continue;

                        }

                        p0Cache.register(variantData, childToParentMap);

                    }

                    pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                    for (short alleleIA = 0; alleleIA < variantInformationA.alleles.length; alleleIA++) {

                        for (short alleleIB = 0; alleleIB < variantInformationB.alleles.length; alleleIB++) {

                            double nA = 0.0;
                            double nB = 0.0;
                            double nAB = 0.0;
                            double n = 0.0;

                            float[] allelePHomA = pHomA[alleleIA];
                            float[] allelePHomB = pHomB[alleleIB];

                            for (int parentI = 0; parentI < allelePHomA.length; parentI++) {

                                float parentAllelePHomA = allelePHomA[parentI];
                                float parentAllelePHomB = allelePHomB[parentI];

                                if (!Float.isNaN(parentAllelePHomA) && !Float.isNaN(parentAllelePHomB)) {

                                    n += 1;

                                    nA += parentAllelePHomA;
                                    nB += parentAllelePHomB;
                                    nAB += parentAllelePHomA * parentAllelePHomB;

                                }
                            }

                            if (nAB * n != nA * nB) {

                                double pAB = nAB / n;
                                double pA = nA / n;
                                double pB = nB / n;

                                double d = pAB - (pA * pB);

                                double r2Value = (d * d) / (pA * (1 - pA) * pB * (1 - pB));

                                writer.writeLine(
                                        variantInformationA.id,
                                        variantInformationA.rsId,
                                        variantInformationA.alleles[alleleIA],
                                        variantInformationB.id,
                                        variantInformationB.rsId,
                                        variantInformationB.alleles[alleleIA],
                                        Double.toString(d),
                                        Double.toString(r2Value)
                                );

                            }
                        }
                    }
                }
            }

            end = Instant.now().getEpochSecond();
            duration = end - start;

            System.out.println("Parsing " + bgenFile + " done (" + duration + " seconds)");

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

    /**
     * Returns a boolean indicating whether the variant has at least two alleles
     * passing the allele frequency.
     *
     * @param variantData The data on the variant to process.
     *
     * @return A boolean indicating whether the variant has at least two alleles
     * passing the allele frequency.
     */
    private static boolean hasAlleles(
            BgenVariantData variantData
    ) {

        int nAlleles = 0;

        for (int allele : variantData.getOrderedAlleles()) {

            double frequency = variantData.getAlleleFrequency(allele);

            if (frequency >= alleleFrequencyThreshold) {

                nAlleles++;

            } else {

                break;

            }
        }

        return nAlleles > 1;

    }

}
