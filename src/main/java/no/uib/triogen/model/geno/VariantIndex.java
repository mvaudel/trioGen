package no.uib.triogen.model.geno;

import java.util.ArrayList;
import java.util.HashMap;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class keeps track of variants giving them a unique index.
 *
 * @author Marc Vaudel
 */
public class VariantIndex {

    /**
     * Variant id to index map.
     */
    private final HashMap<String, Integer> indexMap = new HashMap<>();
    /**
     * List of added identifiers.
     */
    private final ArrayList<String> idList = new ArrayList<>();
    /**
     * Semaphore for the edition of the index.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);
    
    /**
     * Constructor for an empty index.
     */
    public VariantIndex() {
        
        
    }

    /**
     * Constructor taking an array of variants as input.
     *
     * @param variants The array of variants to index.
     */
    public VariantIndex(
            String[] variants
    ) {

        idList.ensureCapacity(variants.length);

        for (int i = 0; i < variants.length; i++) {

            String variantId = variants[i];
            indexMap.put(variantId, i);

            idList.add(variantId);

        }
    }

    /**
     * Adds a variant to the mapping and returns its index.
     *
     * @param variantId The id of the variant.
     *
     * @return The index of the variant.
     */
    public int add(
            String variantId
    ) {

        semaphore.acquire();

        Integer index = indexMap.get(variantId);

        if (index == null) {

            index = idList.size();
            idList.add(variantId);
            indexMap.put(variantId, index);

        }

        semaphore.release();

        return index;

    }

    /**
     * Returns the index for the given variant.
     *
     * @param variantId The id of the variant.
     *
     * @return The corresponding index.
     */
    public int getIndex(
            String variantId
    ) {

        Integer index = indexMap.get(variantId);

        return index != null ? index : add(variantId);

    }

    /**
     * Returns the variants as array indexed according to the mapping.
     *
     * @return The variants as array indexed according to the mapping.
     */
    public String[] getVariants() {

        return idList.toArray(String[]::new);

    }

}
