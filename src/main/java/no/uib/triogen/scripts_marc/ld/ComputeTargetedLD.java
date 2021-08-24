package no.uib.triogen.scripts_marc.ld;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.HashMap;
import java.util.stream.IntStream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.trio_genotypes.VariantList;
import no.uib.triogen.processing.ld.P0Cache;

/**
 * Debugs the LD calculation.
 * 
 * Command: java -Xmx16G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.scripts_marc.ld.ComputeTargetedLD /mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets_ld
 *
 * @author Marc Vaudel
 */
public class ComputeTargetedLD {

    /**
     * The allele frequency threshold to use.
     */
    private static final double alleleFrequencyThreshold = 0.001;

    private static final int maxDistance = 500000;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {
            
            if (args.length != 1) {
                
                throw new IllegalArgumentException("One argument expected.");
                
            }

            File targetsFile = new File(args[0]);
            
            if (!targetsFile.exists()) {
                
                throw new IllegalArgumentException("Targets file '" + args[0] + "' not found.");
                
            }
            
            File trioFile = new File("/mnt/work/marc/moba/run/triogen/pheno/trio");

            System.out.println("Parsing " + trioFile.getAbsolutePath() + ".");

            long start = Instant.now().getEpochSecond();

            ChildToParentMap childToParentMap = ChildToParentMap.fromFile(trioFile);

            long end = Instant.now().getEpochSecond();
            long duration = end - start;

            System.out.println("Parsing " + trioFile + " done (" + duration + " seconds)");

            IntStream.rangeClosed(1, 22)
                    .parallel()
                    .forEach(
                            chromosome -> extractLD(Integer.toString(chromosome), targetsFile, childToParentMap)
                    );

        } catch (Throwable t) {
            t.printStackTrace();
        }

    }

    public static void extractLD(
            String chromosome,
            File targetsFile,
            ChildToParentMap childToParentMap
    ) {

        try {

            VariantList variantList = VariantList.getVariantList(targetsFile, chromosome);

            if (variantList.variantId.length > 0) {

                File bgenFile = new File("/mnt/archive2/marc/moba/phased_bgen/" + chromosome + ".phased.bgen");

                HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(chromosome);
                int defaultMotherPloidy = InheritanceUtils.getDefaultMotherPloidy(chromosome);
                int defaultFatherPloidy = InheritanceUtils.getDefaultFatherPloidy(chromosome);

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

                for (int variantI = 0; variantI < variantList.variantId.length; variantI++) {

                    String rsid = variantList.variantId[variantI];
                    extractLD(chromosome, rsid, bgenIndex, bgenFileReader, childToParentMap);

                }
            }

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

    public static void extractLD(
            String chromosome,
            String rsid,
            BgenIndex bgenIndex,
            BgenFileReader bgenFileReader,
            ChildToParentMap childToParentMap
    ) throws IOException {

        File output = new File("/mnt/work/marc/moba/mobaRun/docs/ld", rsid + "_ld.gz");

        P0Cache p0Cache = new P0Cache(1);

        HashMap<String, double[]> afCache = new HashMap<>();

        System.out.println("Getting LD for " + rsid + ".");

        long start = Instant.now().getEpochSecond();

        ZstdDecompressor decompressor = new ZstdDecompressor();

        int indexA = -1;
        VariantInformation variantInformationA = null;
        double[] af = null;

        for (int i = 0; i < bgenIndex.variantIdArray.length; i++) {

            variantInformationA = bgenIndex.variantInformationArray[i];

            if (variantInformationA.rsid.equals(rsid)) {

                indexA = i;

                break;

            }
        }

        if (indexA == -1) {

            System.out.println(rsid + " not found.");
            return;

        }

        float[][] pHomA = p0Cache.getPHomozygous(variantInformationA.id);

        if (pHomA == null) {

            BgenVariantData variantData = bgenFileReader.getVariantData(indexA);
            variantData.parse(
                    childToParentMap,
                    decompressor
            );

            if (!hasAlleles(variantData)) {

                try (SimpleFileWriter writer = new SimpleFileWriter(output, true)) {

                    writer.writeLine(rsid + " not found.");

                }

                return;

            }

            p0Cache.register(variantData, childToParentMap);

            afCache.put(variantInformationA.id, variantData.getAlleleFrequency());

        }

        pHomA = p0Cache.getPHomozygous(variantInformationA.id);
        int[] allelesA = p0Cache.getOrderedAlleles(variantInformationA.id);
        double[] afA = afCache.get(variantInformationA.id);

        try (SimpleFileWriter writer = new SimpleFileWriter(output, true)) {

            writer.writeLine(
                    "variant_A",
                    "rsid_A",
                    "allele_A",
                    "allele_frequency_A",
                    "variant_B",
                    "rsid_B",
                    "allele_B",
                    "allele_frequency_B",
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

                    afCache.put(variantInformationB.id, variantData.getAlleleFrequency());

                }

                pHomB = p0Cache.getPHomozygous(variantInformationB.id);
                int[] allelesB = p0Cache.getOrderedAlleles(variantInformationB.id);
                double[] afB = afCache.get(variantInformationB.id);

                for (int iA = 0; iA < variantInformationA.alleles.length - 1; iA++) {

                    for (short iB = 0; iB < variantInformationB.alleles.length - 1; iB++) {

                        double nA = 0.0;
                        double nB = 0.0;
                        double nAB = 0.0;
                        double n = 0.0;

                        float[] allelePHomA = pHomA[iA];
                        float[] allelePHomB = pHomB[iB];

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

                            double alleleFrequencyA = afA[allelesA[iA]];
                            double alleleFrequencyB = afB[allelesB[iB]];

                            writer.writeLine(variantInformationA.id,
                                    variantInformationA.rsid,
                                    variantInformationA.alleles[allelesA[iA]],
                                    Double.toString(alleleFrequencyA),
                                    variantInformationB.id,
                                    variantInformationB.rsid,
                                    variantInformationB.alleles[allelesB[iB]],
                                    Double.toString(alleleFrequencyB),
                                    Double.toString(r2Value)
                            );

                        } else {

                            if (variantInformationB.rsid.equals(rsid)) {

                                System.out.println("excluded");

                            }
                        }
                    }
                }
            }
        }

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        System.out.println("Getting LD for " + rsid + " done (" + duration + " seconds)");

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
