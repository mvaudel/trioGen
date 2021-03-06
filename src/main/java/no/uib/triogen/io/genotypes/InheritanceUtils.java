package no.uib.triogen.io.genotypes;

import java.util.HashMap;

/**
 * Utilities for the
 *
 * @author Marc Vaudel
 */
public class InheritanceUtils {

    /**
     * The key for mothers.
     */
    public static final char MOTHER = 'M';
    /**
     * The key for fathers.
     */
    public static final char FATHER = 'F';

    /**
     * The default chromosome inheritance map for children.
     */
    public static final HashMap<String, HashMap<Integer, char[]>> DEFAULT_CHROMOSOME_INHERITANCE_MAP = getDefaultChromosomeInheritanceMap();

    /**
     * Returns the default chromosome inheritance map for children.
     *
     * @return The default chromosome inheritance map for children.
     */
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
        x.put(2, motherFather);

        for (int i = 1; i <= 22; i++) {

            defaultMap.put(Integer.toString(i), simpleDiploid);

        }

        defaultMap.put("X", x);
        defaultMap.put("23", x);

        defaultMap.put("Y", y);
        defaultMap.put("24", y);

        return defaultMap;

    }

    /**
     * Swap mother and father.
     *
     * @param parent The parent to swap.
     *
     * @return The swapped parent.
     */
    public static char swap(
            char parent
    ) {
        switch (parent) {
            case MOTHER:
                return FATHER;
            case FATHER:
                return MOTHER;
            default:
                throw new IllegalArgumentException("Parent " + parent + " not recognized, should be either " + MOTHER + " or " + FATHER + ".");

        }
    }

    /**
     * Returns the default inheritance map for the given contig name.
     *
     * @param contigName The contig name.
     *
     * @return The default inheritance map for the given contig name.
     */
    public static HashMap<Integer, char[]> getDefaultInheritanceMap(
            String contigName
    ) {

        return DEFAULT_CHROMOSOME_INHERITANCE_MAP.get(contigName);

    }

    /**
     * Returns the default maternal ploidy for the given contig.
     * 
     * @param contigName The contig name.
     * 
     * @return The  default maternal ploidy for the given contig.
     */
    public static int getDefaultMotherPloidy(
            String contigName
    ) {

        if (contigName.equals("X")) {

            return 2;

        } else if (contigName.equals("Y")) {

            return 0;

        }

        try {

            int chromosomeNumber = Integer.parseInt(contigName);

            if (chromosomeNumber < 23) {

                return 2;

            } else if (chromosomeNumber == 23) {

                return 2;

            } else if (chromosomeNumber == 24) {

                return 0;

            }

        } catch (Exception e) {

            throw new IllegalArgumentException("Default maternal ploidy for " + contigName + " not implemented.");

        }

        throw new IllegalArgumentException("Default maternal ploidy for " + contigName + " not implemented.");

    }

    /**
     * Returns the default paternal ploidy for the given contig.
     * 
     * @param contigName The contig name.
     * 
     * @return The  default paternal ploidy for the given contig.
     */
    public static int getDefaultFatherPloidy(
            String contigName
    ) {

        if (contigName.equals("X")) {

            return 1;

        } else if (contigName.equals("Y")) {

            return 1;

        }

        try {

            int chromosomeNumber = Integer.parseInt(contigName);

            if (chromosomeNumber < 23) {

                return 1;

            } else if (chromosomeNumber == 23) {

                return 1;

            } else if (chromosomeNumber == 24) {

                return 1;

            }

        } catch (Exception e) {

            throw new IllegalArgumentException("Default paternal ploidy for " + contigName + " not implemented.");

        }

        throw new IllegalArgumentException("Default paternal ploidy for " + contigName + " not implemented.");

    }
}
