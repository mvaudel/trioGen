package no.uib.triogen.model.phenotypes;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.model.covariates.CovariatesHandler;

/**
 * This class parses phenotypes from a pheno file.
 *
 * @author Marc Vaudel
 */
public class PhenotypesHandler {

    /**
     * Name of the column containing the child ids.
     */
    public static String childIdColumn = "child_SentrixID";

    /**
     * Phenotype name to phenotype value map.
     */
    public final ConcurrentHashMap<String, double[]> phenoMap;

    /**
     * Phenotype name to number of valid values. A valid value here is not NA
     * and not infinite.
     */
    public final ConcurrentHashMap<String, Integer> nValidValuesMap;

    /**
     * Mean of phenotypic values.
     */
    public final ConcurrentHashMap<String, Double> phenoMeanMap;

    /**
     * The number of children for which a phenotype was found.
     */
    public final int nChildren;

    /**
     * Constructor from a phenotypes file.
     *
     * @param phenoFile The file containing the phenotypes.
     * @param childrenIds The ids of the children.
     * @param phenoNames The names of the phenotypes to parse, including
     * covariates.
     */
    public PhenotypesHandler(
            File phenoFile,
            String[] childrenIds,
            HashSet<String> phenoNames
    ) {

        phenoMap = new ConcurrentHashMap<>(phenoNames.size());
        nValidValuesMap = new ConcurrentHashMap<>(phenoNames.size());
        phenoMeanMap = new ConcurrentHashMap<>(phenoNames.size());

        phenoNames.parallelStream()
                .forEach(
                        phenoName -> {
                            double[] phenoValues = new double[childrenIds.length];
                            Arrays.fill(phenoValues, Double.NaN);
                            phenoMap.put(phenoName, phenoValues);
                            nValidValuesMap.put(phenoName, 0);
                        }
                );

        HashMap<String, Integer> childIndexMap = new HashMap<>(childrenIds.length);
        int childIndex = 0;

        for (String childId : childrenIds) {

            childIndexMap.put(childId, childIndex++);

        }

        boolean found = false;
        int lineNumber = 1;

        try ( SimpleFileReader phenoReader = SimpleFileReader.getFileReader(phenoFile)) {

            String line = phenoReader.readLine();

            if (line == null) {

                throw new IllegalArgumentException(
                        "Pheno file is empty."
                );

            }

            String[] lineContent = line.split(IoUtils.SEPARATOR);

            int nColumns = lineContent.length;

            if (nColumns < 2) {

                throw new IllegalArgumentException(
                        "Only one column found in phenotypes file. Please verify that the phenotypes file is correctly formatted."
                );

            }

            int idColumnIndex = -1;
            HashMap<String, Integer> columnIndexMap = new HashMap<>(phenoNames.size());

            for (int i = 0; i < lineContent.length; i++) {

                String cellContent = lineContent[i];

                if (cellContent.equals(PhenotypesHandler.childIdColumn)) {

                    idColumnIndex = i;

                } else if (phenoNames.contains(cellContent)) {

                    columnIndexMap.put(cellContent, i);

                }
            }

            if (idColumnIndex == -1) {

                throw new IllegalArgumentException(
                        "The phenotypes file does not contain a \"" + PhenotypesHandler.childIdColumn + "\" column."
                );

            }

            String missingPhenos = phenoNames.stream()
                    .filter(
                            phenoName -> !columnIndexMap.containsKey(phenoName)
                    )
                    .collect(
                            Collectors.joining(", ")
                    );

            if (missingPhenos.length() > 0) {

                throw new IllegalArgumentException(
                        "The following phenotypes were not found in the phenotypes file: " + missingPhenos + "."
                );

            }

            while ((line = phenoReader.readLine()) != null) {

                lineNumber++;

                lineContent = line.split(IoUtils.SEPARATOR);

                if (lineContent.length != nColumns) {

                    throw new IllegalArgumentException(
                            "Unexpected number of columns at line " + lineNumber + " (Expected: " + nColumns + ", Found: " + lineContent.length + ")."
                    );

                }

                String childId = lineContent[idColumnIndex];

                if (childIndexMap.containsKey(childId)) {

                    found = true;

                    childIndex = childIndexMap.get(childId);

                    for (Entry<String, Integer> phenoColumn : columnIndexMap.entrySet()) {

                        String phenoName = phenoColumn.getKey();
                        String valueString = lineContent[phenoColumn.getValue()];

                        if (!valueString.equals("NA") && !valueString.equals("Inf") && !valueString.equals("-Inf")) {

                            double newValue;
                            try {

                                newValue = Double.parseDouble(valueString);

                            } catch (Exception e) {

                                throw new IllegalArgumentException(
                                        "The value for phenotype " + phenoColumn.getKey() + " at for child id " + childId + " (" + valueString + ") could not be parsed as a number."
                                );
                            }

                            double[] phenoValues = phenoMap.get(phenoName);
                            phenoValues[childIndex] = newValue;

                            nValidValuesMap.put(phenoName, nValidValuesMap.get(phenoName) + 1);

                        }
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

    /**
     * Projects phenotypes for the given pheno names using the covariates
     * handler. Removes phenotypes that are not in the given set from the
     * mapping.
     *
     * @param phenoNames
     * @param covariatesHandler
     */
    public void adjustForCovariates(
            String[] phenoNames,
            CovariatesHandler covariatesHandler
    ) {
        
        phenoMap.clear();
        nValidValuesMap.clear();

        Arrays.stream(phenoNames)
                .parallel()
                .forEach(phenoName -> {

                            double[] newValues = covariatesHandler.getProjectedValues(phenoName);
                            
                            double phenoMean = Arrays.stream(newValues)
                                    .sum() / newValues.length;

                            phenoMap.put(phenoName, newValues);
                            nValidValuesMap.put(phenoName, newValues.length);
                            phenoMeanMap.put(phenoName, phenoMean);

                        }
                );
    }
    
    /**
     * Checks that every phenotype has a value and no NaN or infinite value. Throws an exception otherwise.
     */
    public void sanityCheck() {

        phenoMap.entrySet().parallelStream()
                .forEach(
                        entry -> {
                            
                            if (entry.getValue().length == 0) {
                                
                                throw new IllegalArgumentException("Phenotpype " + entry.getKey() + " has no value after filtering and adjustement for covariates.");
                                 
                            }
                            
                            for (double value : entry.getValue()) {
                                
                                if (Double.isNaN(value)) {
                                    
                                    throw new IllegalArgumentException("Phenotpype " + entry.getKey() + " contains NaN value after filtering and adjustement for covariates.");
                                    
                                }
                                if (Double.isInfinite(value)) {
                                    
                                    throw new IllegalArgumentException("Phenotpype " + entry.getKey() + " contains infinite value after filtering and adjustement for covariates.");
                                    
                                }
                            }
                        }
                );
    }
}
