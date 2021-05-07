package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Cache for p0 values.
 *
 * @author Marc Vaudel
 */
public class P0Cache {

    /**
     * The current check-out position for each thread.
     */
    private final int[] checkOutPosition;

    /**
     * Map of the position of each variant.
     */
    private final TreeMap<Integer, ArrayList<String>> positionMap = new TreeMap<>();
    /**
     * Semaphore for the edition of the map.
     */
    private final SimpleSemaphore positionMapSemaphore = new SimpleSemaphore(1);

    /**
     * Map of the probability of being homozygous for each variant.
     */
    private final ConcurrentHashMap<String, float[][]> pHomozygous = new ConcurrentHashMap<>();

    /**
     * Constructor.
     *
     * @param nThreads the number of threads using the cache.
     */
    public P0Cache(
            int nThreads
    ) {

        this.checkOutPosition = new int[nThreads];

    }

    /**
     * Returns the probability of being homozygous for the given variant.
     *
     * @param id The id of the variant.
     *
     * @return A map of the probability of being homozygous per allele and per
     * sample.
     */
    public float[][] getPHomozygous(
            String id
    ) {

        return pHomozygous.get(id);

    }

    /**
     * Releases the given position for the given thread. Once all threads are
     * done at a given position, previous information is cleared from the cache.
     *
     * @param thread The thread number.
     * @param pos The position on the chromosome.
     */
    public void release(int thread, int pos) {

        checkOutPosition[thread] = pos;

        for (int threadPosition : checkOutPosition) {

            if (threadPosition < pos) {

                return;

            }
        }

        positionMapSemaphore.acquire();

        for (Entry<Integer, ArrayList<String>> entry : positionMap.entrySet()) {

            int position = entry.getKey();

            if (position > pos) {

                positionMapSemaphore.release();

                return;

            }

            for (String id : entry.getValue()) {

                pHomozygous.remove(id);

            }
        }

        positionMapSemaphore.release();

    }

    /**
     * Saves the information needed for LD calculation in cache.
     *
     * @param variantData The genotyping data on this variant.
     * @param childToParentMap The child to parent map.
     */
    public void register(
            BgenVariantData variantData,
            ChildToParentMap childToParentMap
    ) {

        VariantInformation variantInformation = variantData.getVariantInformation();

        float[][] variantPHomozygous = new float[variantInformation.alleles.length][2 * childToParentMap.children.length];

        for (int alleleI = 0; alleleI < variantInformation.alleles.length; alleleI++) {

            for (int i = 0; i < childToParentMap.children.length; i++) {

                String childId = childToParentMap.children[i];

                String motherId = childToParentMap.getMother(childId);

                if (variantData.contains(motherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(motherId); z++) {

                        homozygous *= variantData.getProbability(motherId, z, alleleI);

                    }

                    variantPHomozygous[alleleI][i] = homozygous;

                } else {

                    variantPHomozygous[alleleI][i] = Float.NaN;

                }

                String fatherId = childToParentMap.getFather(childId);

                if (variantData.contains(fatherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(fatherId); z++) {

                        homozygous *= variantData.getProbability(fatherId, z, alleleI);

                    }

                    variantPHomozygous[alleleI][i + childToParentMap.children.length] = homozygous;

                } else {

                    variantPHomozygous[alleleI][i + childToParentMap.children.length] = Float.NaN;

                }
            }
        }

        pHomozygous.put(variantInformation.id, variantPHomozygous);
        
        savePosition(variantInformation);

    }

    /**
     * Saves the variant information sorted by position.
     * 
     * @param variantInformation 
     */
    private void savePosition(VariantInformation variantInformation) {

        positionMapSemaphore.acquire();

        ArrayList<String> idsAtPosition = positionMap.get(variantInformation.position);

        if (idsAtPosition == null) {

            idsAtPosition = new ArrayList<>(2);
            positionMap.put(variantInformation.position, idsAtPosition);

        }

        idsAtPosition.add(variantInformation.id);

        positionMapSemaphore.release();

    }
}
