package no.uib.triogen.model.geno;

import java.io.File;
import java.util.ArrayList;
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
     * The chromosome name where to look for.
     */
    public final String[] chromosome;
    /**
     * The position where to start looking for.
     */
    public final int[] start;
    /**
     * The position were to stop looking for.
     */
    public final int[] end;

    /**
     * Constructor.
     *
     * @param variantId the ids of the variant to include
     * @param chromosome the chromosome name where to look for
     * @param start the position where to start looking for
     * @param end the position were to stop looking for
     */
    public VariantList(
            String[] variantId,
            String[] chromosome,
            int[] start,
            int[] end
    ) {
        this.variantId = variantId;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    /**
     * Parses a variant list from a file. First lines starting with '#' are
     * ignored.
     *
     * @param variantFile the file
     *
     * @return an instance of variant list
     */
    public static VariantList getVariantList(
            File variantFile
    ) {

        ArrayList<String> variantIdList = new ArrayList<>();
        ArrayList<String> chromosomeList = new ArrayList<>();
        ArrayList<Integer> startList = new ArrayList<>();
        ArrayList<Integer> endList = new ArrayList<>();

        int lineNumber = 0;

        try (SimpleFileReader phenoReader = SimpleFileReader.getFileReader(variantFile)) {

            String line;
            while ((line = phenoReader.readLine()) != null
                    && (line.trim().length() == 0 || line.charAt(0) == '#')) {

                lineNumber++;

            }

            while ((line = phenoReader.readLine()) != null) {
                
                line = line.trim();

                if (line.length() > 0) {

                    lineNumber++;

                    String[] lineSplit = line.split(IoUtils.SEPARATOR);

                    if (lineSplit.length < 4) {

                        throw new IllegalArgumentException(
                                lineSplit.length + " elements found at line " + lineNumber + " where at least 4 expected. Please make sure that the file is tab-separated.\n" + line
                        );
                    }

                    String variantId = lineSplit[0];
                    String chromosome = lineSplit[1];
                    String startString = lineSplit[2];
                    String endString = lineSplit[3];

                    int start;
                    try {

                        start = Integer.parseInt(startString);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("Start position (" + startString + ") could not be parsed as integer for variant " + variantId + " at line " + lineNumber + ".");

                    }

                    int end;
                    try {

                        end = Integer.parseInt(endString);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("End position (" + endString + ") could not be parsed as integer for variant " + variantId + " at line " + lineNumber + ".");

                    }

                    variantIdList.add(variantId);
                    chromosomeList.add(chromosome);
                    startList.add(start);
                    endList.add(end);

                }
            }
        }

        if (variantIdList.isEmpty()) {

            throw new IllegalArgumentException(
                    "No variant found in " + variantFile + ". Please verify that the header is not commented by \"#\""
            );

        }

        return new VariantList(
                variantIdList.toArray(new String[variantIdList.size()]),
                chromosomeList.toArray(new String[chromosomeList.size()]),
                startList.stream()
                        .mapToInt(a -> a)
                        .toArray(),
                endList.stream()
                        .mapToInt(a -> a)
                        .toArray()
        );
    }
}
