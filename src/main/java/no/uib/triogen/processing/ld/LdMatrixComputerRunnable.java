package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import no.uib.triogen.io.genotypes.iterators.BufferedGenotypesIterator;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.iterators.VariantIterator;
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
     * The buffer.
     */
    private final BufferedGenotypesIterator iterator;
    /**
     * The child ids of the trios to include.
     */
    private final String[] childIds;
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
     * The logger.
     */
    private final Logger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;

    /**
     * Constructor.
     *
     * @param iterator The variant iterator.
     * @param childIds The child ids of the trios to include.
     * @param childToParentMap The map of trios.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param variantIndex The index to use for the variants.
     * @param logger The logger.
     */
    public LdMatrixComputerRunnable(
            BufferedGenotypesIterator iterator,
            String[] childIds,
            ChildToParentMap childToParentMap,
            int maxDistance,
            VariantIndex variantIndex,
            Logger logger
    ) {

        this.iterator = iterator;
        this.childIds = childIds;
        this.childToParentMap = childToParentMap;
        this.maxDistance = maxDistance;
        this.variantIndex = variantIndex;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider genotypesProviderA;
            while ((genotypesProviderA = iterator.next()) != null) {

                GenotypesProvider[] genotypesProviders = iterator.getGenotypesInRange(
                        genotypesProviderA.getContig(),
                        genotypesProviderA.getBp() - maxDistance,
                        genotypesProviderA.getBp() + maxDistance
                );

                int variantIdA = variantIndex.getIndex(genotypesProviderA.getVariantID());
                ArrayList<Integer> variantIds = new ArrayList<>(genotypesProviders.length);
                ArrayList<Double> r2s = new ArrayList<>(genotypesProviders.length);

                for (GenotypesProvider genotypesProviderB : genotypesProviders) {

                    int variantIdB = variantIndex.getIndex(genotypesProviderB.getVariantID());

                    int nAB = 0;
                    int nA = 0;
                    int nB = 0;
                    int n = 2 * childIds.length;

                    for (String childId : childIds) {

                        int[] hA = genotypesProviderA.getH(childToParentMap, childId);
                        int[] hB = genotypesProviderB.getH(childToParentMap, childId);

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

                    if (nA != 0 && nA != n && nB != 0 && nB != n && nAB * n != nA * nB) {

                        double pAB = nAB / n;
                        double pA = nA / n;
                        double pB = nB / n;

                        double d = pAB - (pA * pB);

                        double r2 = (d * d) / (pA * (1 - pA) * pB * (1 - pB));
                        
                        variantIds.add(variantIdB);
                        r2s.add(r2);

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
}
