package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.bgen.VariantIterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.model.trio_genotypes.VariantIndex;

/**
 * Runnable for the LD matrix writer.
 *
 * @author Marc Vaudel
 */
public class LdMatrixComputerRunnable implements Runnable {

    /**
     * The writer.
     */
    private final LdMatrixWriter writer;
    /**
     * The index of the bgen file to process.
     */
    private final BgenIndex bgenIndex;
    /**
     * The reader for the bgen file to process.
     */
    private final BgenFileReader bgenFileReader;
    /**
     * The iterator.
     */
    private final VariantIterator iteratorA;
    /**
     * Cache for the probability of homozygocity.
     */
    private final P0Cache p0Cache;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The bp distance used to compute the ld in sliding windows. A max distance
     * of 10 bp means a sliding window of 20 bp.
     */
    private final int maxDistance;
    /**
     * Index for the variants.
     */
    private final VariantIndex variantIndex;
    /**
     * Index of the thread of this runnable.
     */
    private final int threadIndex;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;
    /**
     * The minimal ld r2 to report (inclusive).
     */
    private final double minR2;

    /**
     * Constructor.
     *
     * @param writer The writer to use.
     * @param iterator The variant iterator.
     * @param bgenIndex The index of the bgen file.
     * @param bgenFileReader The reader for the bgen file.
     * @param childToParentMap The map of trios.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param minR2 The minimal ld r2 to report (inclusive).
     * @param variantIndex The index to use for the variants.
     * @param p0Cache The cache for the probability of homozygocity.
     * @param threadIndex The index of the thread.
     * @param logger The logger.
     */
    public LdMatrixComputerRunnable(
            LdMatrixWriter writer,
            VariantIterator iterator,
            BgenIndex bgenIndex,
            BgenFileReader bgenFileReader,
            ChildToParentMap childToParentMap,
            int maxDistance,
            double minR2,
            VariantIndex variantIndex,
            P0Cache p0Cache,
            int threadIndex,
            SimpleCliLogger logger
    ) {

        this.writer = writer;
        this.iteratorA = iterator;
        this.bgenIndex = bgenIndex;
        this.bgenFileReader = bgenFileReader;
        this.childToParentMap = childToParentMap;
        this.maxDistance = maxDistance;
        this.minR2 = minR2;
        this.variantIndex = variantIndex;
        this.p0Cache = p0Cache;
        this.threadIndex = threadIndex;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            Integer indexA;
            while ((indexA = iteratorA.next()) != null && !canceled) {

                VariantInformation variantInformationA = bgenIndex.variantInformationArray[indexA];

                int variantIdA = variantIndex.getIndex(variantInformationA.id, variantInformationA.rsId);
                    System.out.println("A: " + variantInformationA.id);

                ArrayList<R2> r2s = new ArrayList<>(2);

                float[][] pHomA = p0Cache.getPHomozygous(variantInformationA.id);

                if (pHomA == null) {

                    BgenVariantData variantData = bgenFileReader.getVariantData(indexA);
                    variantData.parse(childToParentMap);

                    p0Cache.register(variantData, childToParentMap);
                    
                }

                pHomA = p0Cache.getPHomozygous(variantInformationA.id);

                VariantIterator iteratorB = new VariantIterator(bgenIndex, variantInformationA.position - maxDistance, variantInformationA.position + maxDistance);

                Integer indexB;
                while ((indexB = iteratorB.next()) != null) {

                    VariantInformation variantInformationB = bgenIndex.variantInformationArray[indexB];
                    int variantIdB = variantIndex.getIndex(variantInformationB.id, variantInformationB.rsId);
//                    System.out.println("B: " + variantInformationA.id);

                    float[][] pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                    if (pHomB == null) {

                        BgenVariantData variantData = bgenFileReader.getVariantData(indexB);
                        variantData.parse(childToParentMap);

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

                                if (r2Value > minR2) {

                                    R2 r2 = new R2(variantIdB, alleleIA, alleleIB, (float) r2Value);
                                    r2s.add(r2);

                                    if (r2Value > 0.8 && variantIdA != variantIdB) {

                                        System.out.println(variantInformationA.rsId + " (" + variantInformationA.id + "@" + alleleIA + ") - " + variantInformationB.rsId + " (" + variantInformationB.id + "@" + alleleIB + "): " + r2Value);

                                        canceled = true;
                                    }

                                }
                            }
                        }
                    }

                    p0Cache.release(threadIndex, variantInformationA.position - maxDistance);

                }
                if (!r2s.isEmpty()) {

                    writer.addVariant(
                            variantIdA,
                            r2s
                    );

                    System.out.println(variantInformationA.rsId + " " + r2s.size());

                } else {

                    System.out.println(variantInformationA.rsId + " No R2");

                }
            }

        } catch (Throwable t) {

            canceled = true;

            logger.logError(
                    Arrays.stream(t.getStackTrace())
                            .map(
                                    element -> element.toString()
                            )
                            .collect(Collectors.joining(" "))
            );

            t.printStackTrace();

        }
    }
}
