package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import no.uib.triogen.io.genotypes.bgen.variant_data.BgenVariantTrioData;
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
     * The current blocking position for each thread.
     */
    private final int[] blockPosition;
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
     * The ordered alleles for each variant.
     */
    private final ConcurrentHashMap<String, int[]> alleles = new ConcurrentHashMap<>();

    /**
     * Constructor.
     *
     * @param nThreads the number of threads using the cache.
     */
    public P0Cache(
            int nThreads
    ) {

        this.blockPosition = new int[nThreads];

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
     * Returns the alleles ordered by allele frequency for the given variant.
     *
     * @param id The id of the variant.
     *
     * @return The alleles orderd by allele freqyency for the given variant.
     */
    public int[] getOrderedAlleles(
            String id
    ) {

        return alleles.get(id);

    }

    /**
     * Releases the given position for the given thread.
     *
     * @param thread The thread number.
     * @param pos The position on the chromosome.
     */
    public void block(
            int thread,
            int pos
    ) {

        blockPosition[thread] = pos;

    }

    /**
     * Releases the given position for the given thread.
     *
     * @param thread The thread number.
     */
    public void release(
            int thread
    ) {

        blockPosition[thread] = -1;

    }

    /**
     * Clear information from the cache.
     *
     * @return Returns the position of the first and last values removed from cache. Null if no value was removed.
     */
    public int[] cleanCache() {
        
        int lastLock = -1;

        for (int threadPosition : blockPosition) {

            if (threadPosition > -1 && 
                    (threadPosition < lastLock || lastLock == -1)) {
                
                        lastLock = threadPosition;

            }
        }
        
        if (lastLock == -1) {
            
            return null;
            
        }

        positionMapSemaphore.acquire();

        HashSet<Integer> positionToRemove = new HashSet<>();

        int firstPosition = -1;

        for (Entry<Integer, ArrayList<String>> entry : positionMap.entrySet()) {

            int position = entry.getKey();

            if (position >= lastLock) {

                break;

            }

            if (firstPosition == -1) {

                firstPosition = position;

            }

            positionToRemove.add(position);

            for (String id : entry.getValue()) {

                pHomozygous.remove(id);
                alleles.remove(id);

            }
        }

        for (int position : positionToRemove) {

            positionMap.remove(position);

        }
        
//        System.out.println("Cache emptied between " + firstPosition + " and " + lastLock);
//        System.out.println("pHomozygous size " + pHomozygous.size());
//        System.out.println("alleles size " + alleles.size());
//        System.out.println("positionMap size " + positionMap.size());

        positionMapSemaphore.release();

        return new int[]{firstPosition, lastLock};

    }

    /**
     * Saves the information needed for LD calculation in cache.
     *
     * @param variantData The genotyping data on this variant.
     * @param childToParentMap The child to parent map.
     */
    public void register(
            BgenVariantTrioData variantData,
            ChildToParentMap childToParentMap
    ) {

        VariantInformation variantInformation = variantData.getVariantInformation();
        int[] orderedAlleles = variantData.getOrderedAlleles();

        float[][] variantPHomozygous = new float[variantInformation.alleles.length - 1][2 * childToParentMap.children.length];

        for (int alleleI = 1; alleleI < orderedAlleles.length; alleleI++) {

            int alleleIndex = orderedAlleles[alleleI];

            for (int childI = 0; childI < childToParentMap.children.length; childI++) {

                String childId = childToParentMap.children[childI];

                String motherId = childToParentMap.getMother(childId);

                if (variantData.contains(motherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(motherId); z++) {

                        homozygous *= variantData.getProbability(motherId, z, alleleIndex);

                    }

                    variantPHomozygous[alleleI - 1][childI] = homozygous;

                } else {

                    variantPHomozygous[alleleI - 1][childI] = Float.NaN;

                }

                String fatherId = childToParentMap.getFather(childId);

                if (variantData.contains(fatherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(fatherId); z++) {

                        homozygous *= variantData.getProbability(fatherId, z, alleleIndex);

                    }

                    variantPHomozygous[alleleI - 1][childI + childToParentMap.children.length] = homozygous;

                } else {

                    variantPHomozygous[alleleI - 1][childI + childToParentMap.children.length] = Float.NaN;

                }
            }
        }

        pHomozygous.put(variantInformation.id, variantPHomozygous);
        alleles.put(variantInformation.id, orderedAlleles);

        savePosition(variantInformation);

    }

    /**
     * Saves the variant information sorted by position.
     *
     * @param variantInformation
     */
    private void savePosition(
            VariantInformation variantInformation
    ) {

        positionMapSemaphore.acquire();

        ArrayList<String> idsAtPosition = positionMap.get(variantInformation.position);

        if (idsAtPosition == null) {

            idsAtPosition = new ArrayList<>(1);
            positionMap.put(variantInformation.position, idsAtPosition);

        }

        idsAtPosition.add(variantInformation.id);

        positionMapSemaphore.release();

    }
}
