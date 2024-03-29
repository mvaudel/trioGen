package no.uib.triogen.model.family;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This class contains mapping between trios.
 *
 * @author Marc Vaudel
 */
public class ChildToParentMap {

    /**
     * Child to mother map.
     */
    private final HashMap<String, String> childToMotherMap;
    /**
     * Child to father map.
     */
    private final HashMap<String, String> childToFatherMap;
    /**
     * Ordered array of children.
     */
    public final String[] children;
    /**
     * Set of all sample ids (child + mother + father).
     */
    public final HashSet<String> sampleIds;

    /**
     * Constructor.
     *
     * @param children an ordered array of all children
     * @param childToMotherMap a child to mother map
     * @param childToFatherMap a child to father map
     */
    public ChildToParentMap(
            String[] children, 
            HashMap<String, String> childToMotherMap, 
            HashMap<String, String> childToFatherMap
    ) {

        this.children = children;
        this.childToFatherMap = childToFatherMap;
        this.childToMotherMap = childToMotherMap;
        
        sampleIds = Arrays.stream(children)
                .collect(
                        Collectors.toCollection(
                                HashSet::new
                        )
                );
        
        sampleIds.addAll(childToFatherMap.values());
        sampleIds.addAll(childToMotherMap.values());
        
    }

    /**
     * Parses the trio map from a file. The file should contain three
     * space-separated columns with one line header. first column child id,
     * second father id, third mother id. NA for missing. Very minimal sanity
     * check is conducted.
     *
     * @param trioFile the file containing the trio map
     * 
     * @return a new instance of a map
     */
    public static ChildToParentMap fromFile(
            File trioFile
    ) {

        TreeSet<String> childIdsSet = new TreeSet<>();
        HashMap<String, String> childToMotherMap = new HashMap<>();
        HashMap<String, String> childToFatherMap = new HashMap<>();

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(trioFile)) {

            String line = reader.readLine();
            String[] lineSplit = line.split(IoUtils.SEPARATOR);

            int lineNumber = 1;

            if (lineSplit.length != 3) {

                throw new IllegalArgumentException("3 elements expected, " + lineSplit.length + " found in line " + lineNumber + " of file " + trioFile + ".\n"
                        + "Plase make sure that the file contains the identifiers of the children, mother, and father, in that order. Columns should be tab-separated.\n" + line);

            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(IoUtils.SEPARATOR);

                if (lineSplit.length != 3) {

                    throw new IllegalArgumentException("3 elements expected, " + lineSplit.length + " found in line " + lineNumber + " of file " + trioFile + ".\n" + line);

                }

                String childId = lineSplit[0];

                if (!childId.equals("NA")) {

                    String motherId = lineSplit[1];
                    String fatherId = lineSplit[2];
                    
                    childIdsSet.add(childId);
                    
                    if (!fatherId.equals("NA")) {
                        
                        childToFatherMap.put(childId, fatherId);
                        
                    }
                    
                    if (!motherId.equals("NA")) {
                        
                        childToMotherMap.put(childId, motherId);
                        
                    }
                }
            }
        }
        
        String[] childIds = childIdsSet.toArray(new String[childIdsSet.size()]);
        
        return new ChildToParentMap(
                childIds, 
                childToMotherMap, 
                childToFatherMap
        );
    }

    /**
     * Returns the father id corresponding to a child.
     *
     * @param child the child id
     *
     * @return the father id
     */
    public String getFather(
            String child
    ) {

        return childToFatherMap.get(child);

    }

    /**
     * Returns the mother id corresponding to a child.
     *
     * @param child the child id
     *
     * @return the mother id
     */
    public String getMother(
            String child
    ) {

        return childToMotherMap.get(child);

    }

}
