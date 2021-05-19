package no.uib.triogen.processing.ld;

import io.airlift.compress.zstd.ZstdCompressor;
import io.airlift.compress.zstd.ZstdDecompressor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Semaphore;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
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
     * The compressor to use.
     */
    private final ZstdCompressor compressor = new ZstdCompressor();
    /**
     * The decompressor to use.
     */
    private final ZstdDecompressor decompressor = new ZstdDecompressor();
    /**
     * The allele frequency threshold to use.
     */
    private final double alleleFrequencyThreshold;
    /**
     * Cache for the variants that have less than two alleles passing the allele
     * frequency threshold.
     */
    private static final HashSet<Integer> excludedVariants = new HashSet<>();
    /**
     * The interval at which this thread should check the cache in number of
     * variants.
     */
    public static final int N_CACHE_FREQ = 10000;
    /**
     * Map of semaphore to synchronize the parsing between threads.
     */
    private final static ConcurrentHashMap<Integer, Semaphore> parseSemaphores = new ConcurrentHashMap<>();

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
     * @param alleleFrequencyThreshold The allele frequency threshold to use.
     * Only variants having at least two alleles passing the threshold will be
     * considered.
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
            double alleleFrequencyThreshold,
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
        this.alleleFrequencyThreshold = alleleFrequencyThreshold;
        this.variantIndex = variantIndex;
        this.p0Cache = p0Cache;
        this.threadIndex = threadIndex;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            int cacheCounter = (new Random()).nextInt(N_CACHE_FREQ);

            Integer indexA;
            while ((indexA = iteratorA.next()) != null && !canceled) {

                if (!excludedVariants.contains(indexA)) {

                    VariantInformation variantInformationA = bgenIndex.variantInformationArray[indexA];

                    float[][] pHomA = p0Cache.getPHomozygous(variantInformationA.id);

                    if (pHomA == null) {

                        Semaphore semaphore = parseSemaphores.get(indexA);

                        if (semaphore != null) {

                            semaphore.acquire();

                            if (excludedVariants.contains(indexA)) {

                                semaphore.release();

                                continue;

                            }

                            pHomA = p0Cache.getPHomozygous(variantInformationA.id);

                        } else {

                            semaphore = new Semaphore(1);
                            semaphore.acquire();
                            parseSemaphores.put(indexA, semaphore);

                        }

                        if (pHomA == null) {

                            BgenVariantData variantData = bgenFileReader.getVariantData(indexA);
                            variantData.parse(
                                    childToParentMap,
                                    decompressor
                            );

                            if (!hasAlleles(variantData)) {

                                if (excludedVariants.size() > 1000000) {

                                    excludedVariants.clear();

                                }

                                excludedVariants.add(indexA);

                                semaphore.release();

                                continue;

                            }

                            p0Cache.register(variantData, childToParentMap);

                            pHomA = p0Cache.getPHomozygous(variantInformationA.id);

                        }

                        semaphore.release();

                    }

                    int variantIdA = variantIndex.getIndex(variantInformationA.id, variantInformationA.rsId);

                    ArrayList<R2> r2s = new ArrayList<>(2);

                    int[] allelesA = p0Cache.getOrderedAlleles(variantInformationA.id);

                    VariantIterator iteratorB = new VariantIterator(bgenIndex, variantInformationA.position - maxDistance, variantInformationA.position + maxDistance);

                    Integer indexB;
                    while ((indexB = iteratorB.next()) != null) {

                        if (!excludedVariants.contains(indexB)) {

                            VariantInformation variantInformationB = bgenIndex.variantInformationArray[indexB];
                            int variantIdB = variantIndex.getIndex(variantInformationB.id, variantInformationB.rsId);

                            float[][] pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                            if (pHomB == null) {

                                Semaphore semaphore = parseSemaphores.get(indexB);

                                if (semaphore != null) {

                                    semaphore.acquire();

                                    if (excludedVariants.contains(indexB)) {

                                        semaphore.release();

                                        continue;

                                    }

                                    pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                                } else {

                                    semaphore = new Semaphore(1);
                                    semaphore.acquire();
                                    parseSemaphores.put(indexB, semaphore);

                                }

                                if (pHomB == null) {

                                    BgenVariantData variantData = bgenFileReader.getVariantData(indexB);
                                    variantData.parse(
                                            childToParentMap,
                                            decompressor
                                    );

                                    if (!hasAlleles(variantData)) {

                                        excludedVariants.add(indexB);

                                        semaphore.release();

                                        continue;

                                    }

                                    p0Cache.register(variantData, childToParentMap);

                                    pHomB = p0Cache.getPHomozygous(variantInformationB.id);

                                    semaphore.release();

                                }
                            }

                            int[] allelesB = p0Cache.getOrderedAlleles(variantInformationB.id);

                            for (int alleleIA = 0; alleleIA < allelesA.length - 1; alleleIA++) {

                                for (int alleleIB = 0; alleleIB < allelesB.length - 1; alleleIB++) {

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

                                            R2 r2 = new R2(variantIdB, (short) allelesA[alleleIA + 1], (short) allelesB[alleleIB + 1], (float) r2Value);
                                            r2s.add(r2);

                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (!r2s.isEmpty()) {

                        writer.addVariant(
                                variantIdA,
                                r2s,
                                compressor
                        );
                    }

                    if (++cacheCounter >= N_CACHE_FREQ) {

                        p0Cache.releaseAndEmptyCache(threadIndex, variantInformationA.position - maxDistance);

                        cacheCounter = 0;

                    } else {

                        p0Cache.release(threadIndex, variantInformationA.position - maxDistance);

                    }
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

    /**
     * Returns a boolean indicating whether the variant has at least two alleles
     * passing the allele frequency.
     *
     * @param variantData The data on the variant to process.
     *
     * @return A boolean indicating whether the variant has at least two alleles
     * passing the allele frequency.
     */
    private boolean hasAlleles(
            BgenVariantData variantData
    ) {

        int nAlleles = 0;
        int[] orderedAlleles = variantData.getOrderedAlleles();

        for (int i = 1; i < orderedAlleles.length; i++) {

            double frequency = variantData.getAlleleFrequency(orderedAlleles[i]);

            if (frequency >= alleleFrequencyThreshold) {

                nAlleles++;

            } else {

                break;

            }
        }

        return nAlleles > 0;

    }
}
