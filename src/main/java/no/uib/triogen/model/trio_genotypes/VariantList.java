package no.uib.triogen.model.trio_genotypes;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This class provides variants to look for.
 *
 * @author Marc Vaudel
 */
public class VariantList {

    /**
     * The ids of the variant to include.
     */
    public final String[] variantId;
    /**
     * The name of the contig.
     */
    public final String[] chromosome;
    /**
     * The position where to look for.
     */
    public final int[] position;
    /**
     * Set of the ids of the variants in the list.
     */
    private final HashMap<String, Integer> variantIdsMap;
    /**
     * Map of the variants indexed by contig and position.
     */
    private HashMap<String, HashMap<Integer, ArrayList<Integer>>> positionToVariantMap;
    /**
     * The distance used to make the index.
     */
    private int indexDistance;

    /**
     * Constructor.
     *
     * @param variantId The ids of the variant to include.
     * @param contig The name of the contig.
     * @param position The position of the variant.
     */
    public VariantList(
            String[] variantId,
            String[] contig,
            int[] position
    ) {
        this.variantId = variantId;
        this.chromosome = contig;
        this.position = position;

        variantIdsMap = new HashMap<>(variantId.length);

        for (int i = 0; i < variantId.length; i++) {

            variantIdsMap.put(variantId[i], i);

        }
    }

    /**
     * Parses a variant list from a file. First lines starting with '#' are
     * ignored.
     *
     * @param variantFile The file containing the variant information.
     *
     * @return an instance of variant list
     */
    public static VariantList getVariantList(
            File variantFile
    ) {
        
        return getVariantList(variantFile, null);
        
    }

    /**
     * Parses a variant list from a file. First lines starting with '#' are
     * ignored.
     *
     * @param variantFile The file containing the variant information.
     * @param chromosome The chromosome, ignored if null.
     *
     * @return an instance of variant list
     */
    public static VariantList getVariantList(
            File variantFile,
            String chromosome
    ) {

        ArrayList<String> variantIdList = new ArrayList<>();
        ArrayList<String> contigList = new ArrayList<>();
        ArrayList<Integer> positionList = new ArrayList<>();

        int lineNumber = 0;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(variantFile)) {

            String line;
            while ((line = reader.readLine()) != null
                    && (line.trim().length() == 0 || line.charAt(0) == '#')) {

                lineNumber++;

            }

            while ((line = reader.readLine()) != null) {

                line = line.trim();

                if (line.length() > 0) {

                    lineNumber++;

                    String[] lineSplit = line.split(IoUtils.SEPARATOR);

                    if (lineSplit.length < 3) {

                        throw new IllegalArgumentException(
                                lineSplit.length + " elements found at line " + lineNumber + " where at least 3 expected. Please make sure that the file is tab-separated.\n" + line
                        );
                    }

                    String variantId = lineSplit[0];
                    String contig = lineSplit[1];

                    if (chromosome == null || contig.equals(chromosome)) {

                        String positionString = lineSplit[2];

                        int position;
                        try {

                            position = Integer.parseInt(positionString);

                        } catch (Exception e) {

                            throw new IllegalArgumentException("Position (" + positionString + ") could not be parsed as integer for variant " + variantId + " at line " + lineNumber + ".");

                        }

                        variantIdList.add(variantId);
                        contigList.add(contig);
                        positionList.add(position);

                    }
                }
            }
        }

        return new VariantList(
                variantIdList.toArray(new String[variantIdList.size()]),
                contigList.toArray(new String[contigList.size()]),
                positionList.stream()
                        .mapToInt(a -> a)
                        .toArray()
        );
    }

    /**
     * Returns a boolean indicating whether the variant id is in the list.
     *
     * @param variantId The id of the variant of interest.
     *
     * @return A boolean indicating whether the variant id is in the list.
     */
    public boolean contains(
            String variantId
    ) {

        return variantIdsMap.containsKey(variantId);

    }

    /**
     * Returns the index of the given variant id in the list.
     *
     * @param variantId The id of the variant of interest.
     *
     * @return The index of the given variant id, null if not found.
     */
    public Integer getIndex(
            String variantId
    ) {

        return variantIdsMap.get(variantId);

    }

    /**
     * Indexes the given variant list using the given reference distance.
     *
     * @param indexDistance The reference distance.
     */
    public void index(
            int indexDistance
    ) {

        if (indexDistance == 0) {

            indexDistance = 1;

        }

        positionToVariantMap = new HashMap<>(4);
        this.indexDistance = indexDistance;

        for (int variantI = 0; variantI < variantId.length; variantI++) {

            int positionI = this.position[variantI];
            String contigI = this.chromosome[variantI];

            HashMap<Integer, ArrayList<Integer>> contigIndex = positionToVariantMap.get(contigI);

            if (contigIndex == null) {

                contigIndex = new HashMap<>(2);
                positionToVariantMap.put(contigI, contigIndex);

            }

            int refIndex = getPositionIndex(positionI, indexDistance);

            ArrayList<Integer> variantsAtIndex = contigIndex.get(refIndex);

            if (variantsAtIndex == null) {

                variantsAtIndex = new ArrayList<>(1);
                contigIndex.put(refIndex, variantsAtIndex);

            }

            variantsAtIndex.add(variantI);

        }
    }

    /**
     * Indicates whether a variant at the given position on the given contig
     * should be included. The distance tolerance used is the one that was used
     * to create the index.
     *
     * @param contig The name of the contig.
     * @param position The position on the contig.
     *
     * @return A boolean indicating whether a variant at the given position on
     * the given contig should be included.
     */
    public boolean include(
            String contig, 
            int position
    ) {

        HashMap<Integer, ArrayList<Integer>> contigIndex = positionToVariantMap.get(contig);

        if (contigIndex != null) {

            int index = getPositionIndex(position, indexDistance);

            for (int testIndex = index - 1; testIndex <= index + 2; testIndex++) {

                ArrayList<Integer> tempIndexes = contigIndex.get(testIndex);

                if (tempIndexes != null) {

                    for (int variantI : tempIndexes) {

                        if (Math.abs(this.position[variantI] - position) <= indexDistance) {

                            return true;

                        }
                    }
                }
            }
        }

        return false;

    }

    /**
     * Returns the index for the given position.
     *
     * @param position The position on the contig.
     * @param indexDistance The reference distance to use in the index.
     *
     * @return The index for the given position.
     */
    public static int getPositionIndex(
            int position,
            int indexDistance
    ) {

        return position / indexDistance;

    }
}
