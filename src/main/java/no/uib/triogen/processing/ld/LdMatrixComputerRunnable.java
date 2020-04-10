package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import no.uib.triogen.io.genotypes.iterators.BufferedGenotypesIterator;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.log.Logger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.VariantIndex;

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
     * The buffer.
     */
    private final BufferedGenotypesIterator iterator;
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
    private final Logger logger;
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
            BufferedGenotypesIterator iterator,
            ChildToParentMap childToParentMap,
            int maxDistance,
            double minR2,
            boolean hardCalls,
            VariantIndex variantIndex,
            boolean testIteration,
            Logger logger
    ) {

        this.writer = writer;
        this.iterator = iterator;
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
            while ((genotypesProviderA = iterator.next()) != null && !canceled) {

                GenotypesProvider[] genotypesProviders = iterator.getGenotypesInRange(
                        genotypesProviderA.getContig(),
                        genotypesProviderA.getBp() - maxDistance,
                        genotypesProviderA.getBp() + maxDistance
                );

                int variantIdA = variantIndex.getIndex(genotypesProviderA.getVariantID());
                ArrayList<Integer> variantIds = new ArrayList<>(genotypesProviders.length);
                ArrayList<Double> r2s = new ArrayList<>(genotypesProviders.length);

                if (!testIteration) {

                    for (GenotypesProvider genotypesProviderB : genotypesProviders) {

                        int variantIdB = variantIndex.getIndex(genotypesProviderB.getVariantID());

                        if (variantIdA != variantIdB) {

                            double nAB = 0.0;
                            double nA = 0.0;
                            double nB = 0.0;
                            double n = 2 * childToParentMap.children.length;

                            for (String childId : childToParentMap.children) {

                                if (!hardCalls) {

                                    String motherId = childToParentMap.getMother(childId);

                                    float[] dosagesA = genotypesProviderA.getDosages(motherId);
                                    float[] dosagesB = genotypesProviderA.getDosages(motherId);

                                    float pA0 = dosagesA[0];
                                    float pB0 = dosagesB[0];

                                    nA += pA0;
                                    nB += pB0;
                                    nAB += pA0 * pB0;

                                    String fatherId = childToParentMap.getFather(childId);

                                    dosagesA = genotypesProviderA.getDosages(fatherId);
                                    dosagesB = genotypesProviderA.getDosages(fatherId);

                                    pA0 = dosagesA[0];
                                    pB0 = dosagesB[0];

                                    nA += pA0;
                                    nB += pB0;
                                    nAB += pA0 * pB0;

                                } else {

                                    short[] hA = genotypesProviderA.getH(childToParentMap, childId);
                                    short[] hB = genotypesProviderB.getH(childToParentMap, childId);

                                    boolean a = hA[0] == 0 && hA[1] == 0;
                                    boolean b = hB[0] == 0 && hB[1] == 0;

                                    if (a) {
                                        nA++;
                                        if (b) {
                                            nAB++;
                                        }
                                    }

                                    if (b) {
                                        nB++;
                                    }

                                    a = hA[2] == 0 && hA[3] == 0;
                                    b = hB[2] == 0 && hB[3] == 0;

                                    if (a) {
                                        nA++;
                                        if (b) {
                                            nAB++;
                                        }
                                    }

                                    if (b) {
                                        nB++;
                                    }
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
                                    
                                    if (genotypesProviderA.getVariantID().equals("rs79847867")) {
                                        
                                        System.out.println("r2:" + r2 + " pA:" + pA + " pB:" + pB + " pAB:" + pAB + " n:" + n);
                                        
                                    }

                                }
                            }
                        }
                    }
                }

                iterator.releaseMinBp(genotypesProviderA.getContig(), genotypesProviderA.getBp());

                if (!variantIds.isEmpty()) {

                    writer.addVariant(
                            variantIdA,
                            variantIds,
                            r2s
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
}
