package no.uib.triogen.model.trio_genotypes;

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
     * List of added rsIds.
     */
    private final ArrayList<String> rsidList = new ArrayList<>();
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
     * Adds a variant to the mapping and returns its index.
     *
     * @param variantId The id of the variant.
     * @param rsId The rsid of the variant.
     *
     * @return The index of the variant.
     */
    public int add(
            String variantId,
            String rsId
    ) {

        semaphore.acquire();

        Integer index = indexMap.get(variantId);

        if (index == null) {

            index = idList.size();
            idList.add(variantId);
            indexMap.put(variantId, index);

            if (rsId != null) {

                rsidList.add(rsId);

            } else {

                rsidList.add("");

            }
        }

        semaphore.release();

        return index;

    }

    /**
     * Returns the index for the given variant, null if not found.
     *
     * @param variantId The id of the variant.
     * @param rsId The rsid of the variant.
     *
     * @return The corresponding index.
     */
    public int getIndex(
            String variantId,
            String rsId
    ) {

        Integer index = indexMap.get(variantId);

        return index != null ? index : add(variantId, rsId);

    }

    /**
     * Returns the ids of the variants as array indexed according to the
     * mapping.
     *
     * @return The ids of the variants as array indexed according to the
     * mapping.
     */
    public String[] getVariantIds() {

        return idList.stream().toArray(String[]::new);

    }

    /**
     * Returns the rsids of the variants as array indexed according to the
     * mapping.
     *
     * @return The rsids of the variants as array indexed according to the
     * mapping.
     */
    public String[] getRsIds() {

        return rsidList.stream().toArray(String[]::new);

    }

}
