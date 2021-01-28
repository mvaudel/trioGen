package no.uib.triogen.io.genotypes;

import java.util.HashMap;

/**
 * Utilities for the
 *
 * @author Marc Vaudel
 */
public class InheritanceUtils {

    public static final char MOTHER = 'M';
    public static final char FATHER = 'F';

    public static final HashMap<String, HashMap<Integer, char[]>> DEFAULT_CHROMOSOME_INHERITANCE_MAP = getDefaultChromosomeInheritanceMap();

    private static HashMap<String, HashMap<Integer, char[]>> getDefaultChromosomeInheritanceMap() {

        HashMap<String, HashMap<Integer, char[]>> defaultMap = new HashMap<>(25);

        char[] motherFather = new char[]{MOTHER, FATHER};
        char[] mother = new char[]{MOTHER};
        char[] father = new char[]{FATHER};

        HashMap<Integer, char[]> simpleDiploid = new HashMap<>(1);
        simpleDiploid.put(2, motherFather);

        HashMap<Integer, char[]> x = new HashMap<>(2);
        x.put(1, mother);
        x.put(2, motherFather);

        HashMap<Integer, char[]> y = new HashMap<>(1);
        y.put(1, father);

        for (int i = 1; i <= 22; i++) {

            defaultMap.put(Integer.toString(i), simpleDiploid);

        }

        defaultMap.put("X", x);
        defaultMap.put("23", x);

        defaultMap.put("Y", y);
        defaultMap.put("24", y);

        return defaultMap;

    }

}
