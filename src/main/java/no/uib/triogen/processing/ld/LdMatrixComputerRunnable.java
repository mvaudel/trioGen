package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.zip.Deflater;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.trio_genotypes.VariantIndex;
import no.uib.triogen.io.genotypes.WindowGenotypesIterator;

/**
 * Runnable for the LD matrix writer.
 *
 * @author Marc Vaudel
 */
public class LdMatrixComputerRunnable implements Runnable, AutoCloseable {

    /**
     * The writer.
     */
    private final LdMatrixWriter writer;
    /**
     * The buffer.
     */
    private final WindowGenotypesIterator iteratorA;
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
     * Boolean indicating whether hard calls should be used.
     */
    private final boolean hardCalls;
    /**
     * Index for the variants.
     */
    private final VariantIndex variantIndex;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;
    /**
     * Boolean indicating whether the runnable should only iterate through the
     * variants and not compute LD calculations.
     */
    private final boolean testIteration;
    /**
     * The minimal ld r2 to report (inclusive).
     */
    private final double minR2;
    /**
     * The deflater to compress parts of the file.
     */
    private final Deflater deflater = new Deflater(Deflater.BEST_COMPRESSION, true);

    /**
     * Constructor.
     *
     * @param writer The writer to use.
     * @param iterator The variant iterator.
     * @param childToParentMap The map of trios.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param minR2 The minimal ld r2 to report (inclusive).
     * @param hardCalls Boolean indicating whether hard calls should be used
     * instead of dosages.
     * @param variantIndex The index to use for the variants.
     * @param testIteration Boolean indicating whether the runnable should only
     * iterate through the variants and not compute LD calculations.
     * @param logger The logger.
     */
    public LdMatrixComputerRunnable(
            LdMatrixWriter writer,
            WindowGenotypesIterator iterator,
            ChildToParentMap childToParentMap,
            int maxDistance,
            double minR2,
            boolean hardCalls,
            VariantIndex variantIndex,
            boolean testIteration,
            SimpleCliLogger logger
    ) {

        this.writer = writer;
        this.iteratorA = iterator;
        this.childToParentMap = childToParentMap;
        this.maxDistance = maxDistance;
        this.minR2 = minR2;
        this.hardCalls = hardCalls;
        this.variantIndex = variantIndex;
        this.testIteration = testIteration;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider genotypesProviderA;
            while ((genotypesProviderA = iteratorA.next()) != null && !canceled) {

                int variantIdA = variantIndex.getIndex(genotypesProviderA.getVariantID());
                ArrayList<Integer> variantIds = new ArrayList<>();
                ArrayList<Double> r2s = new ArrayList<>();

                VariantIterator iteratorB = iteratorA.getGenotypesInRange(
                        genotypesProviderA.getContig(),
                        Math.max(genotypesProviderA.getBp() - maxDistance, 0),
                        genotypesProviderA.getBp() + maxDistance
                );

                if (!testIteration) {

                    GenotypesProvider genotypesProviderB;
                    while ((genotypesProviderB = iteratorB.next()) != null) {

                        int variantIdB = variantIndex.getIndex(genotypesProviderB.getVariantID());

                        if (!hardCalls) {

                            double nA = genotypesProviderA.getParentsDosageP0Cache();
                            double nB = genotypesProviderB.getParentsDosageP0Cache();
                            double n = 2 * childToParentMap.children.length;

                            if (nA > 0 && nA < n || nB > 0 && nB < n) {

                                float[] p0sA = genotypesProviderA.getParentsDosageP0sCache();
                                float[] p0sB = genotypesProviderB.getParentsDosageP0sCache();
                                double nAB = 0.0;

                                for (int i = 0; i < p0sA.length; i++) {

                                    nAB += p0sA[i] * p0sB[i];

                                }

                                if (nAB * n != nA * nB) {

                                    double pAB = nAB / n;
                                    double pA = nA / n;
                                    double pB = nB / n;

                                    double d = pAB - (pA * pB);

                                    double r2 = (d * d) / (pA * (1 - pA) * pB * (1 - pB));

                                    if (r2 > minR2) {

                                        variantIds.add(variantIdB);
                                        r2s.add(r2);

                                    }
                                }
                            }
                        } else {

                            double nA = genotypesProviderA.getParentsGenotypeP0Cache();
                            double nB = genotypesProviderB.getParentsGenotypeP0Cache();
                            double n = 2 * childToParentMap.children.length;

                            if (nA < n && nA > 0 || nB < n && nB > 0) {

                                boolean[] p0sA = genotypesProviderA.getParentsGenotypeP0sCache();
                                boolean[] p0sB = genotypesProviderB.getParentsGenotypeP0sCache();
                                double nAB = 0.0;

                                for (int i = 0; i < p0sA.length; i++) {

                                    if (p0sA[i] && p0sB[i]) {

                                        nAB += 1.0;

                                    }
                                }

                                if (nA >= 0.0 && nA <= n && nB >= 0 && nB <= n && nAB * n != nA * nB) {

                                    double pAB = nAB / n;
                                    double pA = nA / n;
                                    double pB = nB / n;

                                    double d = pAB - (pA * pB);

                                    double r2 = (d * d) / (pA * (1 - pA) * pB * (1 - pB));

                                    if (r2 > minR2) {

                                        variantIds.add(variantIdB);
                                        r2s.add(r2);

                                    }
                                }
                            }
                        }
                    }
                }

                iteratorA.releaseMinBp(genotypesProviderA.getContig(), genotypesProviderA.getBp());

                if (!variantIds.isEmpty()) {

                    writer.addVariant(
                            variantIdA,
                            variantIds,
                            r2s,
                            deflater
                    );
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

    @Override
    public void close() throws Exception {

        deflater.end();

    }
}
