package no.uib.triogen.processing.ld;

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
public class LdMatrixWriterRunnable implements Runnable {

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
     * @param childToParentMap The map of trios.
     * @param maxDistance The maximal number of bp to allow between variants.
     * @param variantIndex The index to use for the variants.
     * @param logger The logger.
     */
    public LdMatrixWriterRunnable(
            BufferedGenotypesIterator iterator,
            ChildToParentMap childToParentMap,
            int maxDistance,
            VariantIndex variantIndex,
            Logger logger
    ) {

        this.iterator = iterator;
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

                for (GenotypesProvider genotypesProviderB : genotypesProviders) {

                    int variantIdB = variantIndex.getIndex(genotypesProviderB.getVariantID());

                    int nAB = 0;
                    int nA = 0;
                    int nB = 0;
                    int n = 2 * childToParentMap.children.length;

                    for (String childId : childToParentMap.children) {

                        int[] hA = genotypesProviderA.getH(childToParentMap, childId);
                        int[] hB = genotypesProviderB.getH(childToParentMap, childId);

                        boolean a = hA[0] == 0 && hA[1] == 0 && hA[2] == 0 && hA[3] == 0;
                        boolean b = hB[0] == 0 && hB[1] == 0 && hB[2] == 0 && hB[3] == 0;

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

                    double pAB = nAB / n;
                    double pA = nA / n;
                    double pB = nB / n;

                    double d = pAB - (pA * pB);

                    double r2 = (d * d) / (pA * (1 - pA) * pB * (1 - pB));

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
