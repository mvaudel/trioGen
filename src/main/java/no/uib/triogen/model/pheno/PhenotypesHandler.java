package no.uib.triogen.model.pheno;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This class parses phenotypes from a pheno file.
 *
 * @author Marc Vaudel
 */
public class PhenotypesHandler {

    /**
     * Name of the column containing the child ids.
     */
    public final static String childIdColumn = "child_SentrixID";

    /**
     * Phenotype name to child id to phenotype value map.
     */
    public final HashMap<String, double[]> phenoMap;

    /**
     * The number of children for which a phenotype was found.
     */
    public final int nChildren;

    /**
     * Constructor from a phenotypes file.
     *
     * @param phenoFile a file containing the phenotypes
     * @param childrenIds childrenIds
     * @param phenoNames the names of the phenotypes to parse
     */
    public PhenotypesHandler(
            File phenoFile,
            TreeSet<String> childrenIds,
            String[] phenoNames
    ) {

        phenoMap = new HashMap<>(phenoNames.length);

        for (String phenoName : phenoNames) {

            double[] phenoValues = new double[childrenIds.size()];
            Arrays.fill(phenoValues, Double.NaN);
            phenoMap.put(phenoName, phenoValues);

        }

        HashMap<String, Integer> childIndexMap = new HashMap<>(childrenIds.size());
        int childIndex = 0;

        for (String childId : childrenIds) {

            childIndexMap.put(childId, childIndex++);

        }

        boolean found = false;
        int lineNumber = 1;

        try (SimpleFileReader phenoReader = SimpleFileReader.getFileReader(phenoFile)) {

            String line = phenoReader.readLine();

            if (line == null) {

                throw new IllegalArgumentException(
                        "Pheno file is empty."
                );

            }

            String[] lineContent = line.split(Utils.separator);

            int nColumns = lineContent.length;

            if (nColumns < 2) {

                throw new IllegalArgumentException(
                        "Only one column found in pheno file. Please verify that the pheno file is correctly formatted."
                );

            }

            HashSet<String> phenoNamesSet = Arrays.stream(phenoNames)
                    .collect(Collectors.toCollection(HashSet::new));

            int idColumnIndex = -1;
            HashMap<String, Integer> phenoColumnIndexMap = new HashMap<>(phenoNames.length);

            for (int i = 0; i < lineContent.length; i++) {

                String cellContent = lineContent[i];

                if (cellContent.equals(PhenotypesHandler.childIdColumn)) {

                    idColumnIndex = i;

                } else if (phenoNamesSet.contains(cellContent)) {

                    phenoColumnIndexMap.put(cellContent, i);

                }
            }

            if (idColumnIndex == -1) {

                throw new IllegalArgumentException(
                        "The phenotypes file does not contain a \"" + PhenotypesHandler.childIdColumn + "\" column."
                );

            }

            String[] missingPhenos = Arrays.stream(phenoNames)
                    .filter(
                            phenoName -> !phenoColumnIndexMap.containsKey(phenoName)
                    )
                    .toArray(String[]::new);

            if (missingPhenos.length > 0) {

                throw new IllegalArgumentException(
                        "The following phenotypes were not found in the phenotypes file: " + String.join(", ", missingPhenos) + "."
                );

            }

            while ((line = phenoReader.readLine()) != null) {

                lineNumber++;

                lineContent = line.split(Utils.separator);

                if (lineContent.length != nColumns) {

                    throw new IllegalArgumentException(
                            "Unexpected number of columns at line " + lineNumber + " (Expected: " + nColumns + ", Found: " + lineContent.length + ")."
                    );

                }

                String childId = lineContent[idColumnIndex];

                if (childIndexMap.containsKey(childId)) {

                    found = true;

                    childIndex = childIndexMap.get(childId);

                    for (Entry<String, Integer> phenoColumn : phenoColumnIndexMap.entrySet()) {

                        double[] phenoValues = phenoMap.get(phenoColumn.getKey());

                        String valueString = lineContent[phenoColumn.getValue()];

                        double newValue;
                        try {

                            newValue = Double.parseDouble(valueString);

                        } catch (Exception e) {

                            throw new IllegalArgumentException(
                                    "The value for phenotype " + phenoColumn.getKey() + " at for child id " + childId + " (" + valueString + ") could not be parsed as a number."
                            );
                        }

                        phenoValues[childIndex] = newValue;

                    }
                }
            }
        }

        nChildren = lineNumber - 1;

        if (!found) {

            throw new IllegalArgumentException(
                    "No child identifier in the pheno file matched a child identifier in the trio file. As a result, no phenotype is available for the analysis. Please check that the identifiers are correct."
            );

        }
    }
}
