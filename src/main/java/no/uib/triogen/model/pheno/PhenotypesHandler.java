package no.uib.triogen.model.pheno;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

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
     * Phenotype name to child id to phenotype value map.
     */
    public final HashMap<String, double[]> phenoMap;

    /**
     * Phenotype name to number of valid values. A valid value here is not NA
     * and not infinite.
     */
    public final HashMap<String, Integer> nValidValuesMap;

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
     * @param covariates the names of the covariates to parse
     */
    public PhenotypesHandler(
            File phenoFile,
            TreeSet<String> childrenIds,
            String[] phenoNames,
            String[] covariates
    ) {

        phenoMap = new HashMap<>(phenoNames.length);
        nValidValuesMap = new HashMap<>(phenoNames.length);

        for (String phenoName : phenoNames) {

            double[] phenoValues = new double[childrenIds.size()];
            Arrays.fill(phenoValues, Double.NaN);
            phenoMap.put(phenoName, phenoValues);

            nValidValuesMap.put(phenoName, 0);

        }

        double[][] covariatesX = new double[childrenIds.size()][covariates.length];

        for (int i = 0; i < childrenIds.size(); i++) {
            for (int j = 0; j < covariates.length; j++) {

                covariatesX[i][j] = Double.NaN;

            }
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

            String[] lineContent = line.split(IoUtils.separator);

            int nColumns = lineContent.length;

            if (nColumns < 2) {

                throw new IllegalArgumentException(
                        "Only one column found in pheno file. Please verify that the pheno file is correctly formatted."
                );

            }

            HashSet<String> phenoNamesSet = Arrays.stream(phenoNames)
                    .collect(Collectors.toCollection(HashSet::new));
            HashSet<String> covariatesSet = Arrays.stream(covariates)
                    .collect(Collectors.toCollection(HashSet::new));

            int idColumnIndex = -1;
            HashMap<String, Integer> phenoColumnIndexMap = new HashMap<>(phenoNames.length);
            TreeMap<String, Integer> covariatesColumnIndexMap = new TreeMap<>();

            for (int i = 0; i < lineContent.length; i++) {

                String cellContent = lineContent[i];

                if (cellContent.equals(PhenotypesHandler.childIdColumn)) {

                    idColumnIndex = i;

                } else if (phenoNamesSet.contains(cellContent)) {

                    phenoColumnIndexMap.put(cellContent, i);

                } else if (covariatesSet.contains(cellContent)) {

                    covariatesColumnIndexMap.put(cellContent, i);

                }
            }

            if (idColumnIndex == -1) {

                throw new IllegalArgumentException(
                        "The phenotypes file does not contain a \"" + PhenotypesHandler.childIdColumn + "\" column."
                );

            }

            String missingPhenos = Arrays.stream(phenoNames)
                    .filter(
                            phenoName -> !phenoColumnIndexMap.containsKey(phenoName)
                    )
                    .collect(
                            Collectors.joining(", ")
                    );

            if (missingPhenos.length() > 0) {

                throw new IllegalArgumentException(
                        "The following phenotypes were not found in the phenotypes file: " + missingPhenos + "."
                );

            }

            String missingCovariates = Arrays.stream(covariates)
                    .filter(
                            phenoName -> !covariatesColumnIndexMap.containsKey(phenoName)
                    )
                    .collect(
                            Collectors.joining(", ")
                    );

            if (missingCovariates.length() > 0) {

                throw new IllegalArgumentException(
                        "The following covariates were not found in the phenotypes file: " + missingCovariates + "."
                );

            }

            while ((line = phenoReader.readLine()) != null) {

                lineNumber++;

                lineContent = line.split(IoUtils.separator);

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

                    int i = 0;
                    for (Entry<String, Integer> covariateColumn : covariatesColumnIndexMap.entrySet()) {

                        String valueString = lineContent[covariateColumn.getValue()];

                        if (!valueString.equals("NA") && !valueString.equals("Inf") && !valueString.equals("-Inf")) {

                            double newValue;
                            try {

                                newValue = Double.parseDouble(valueString);

                            } catch (Exception e) {

                                throw new IllegalArgumentException(
                                        "The value for covariate " + covariateColumn.getKey() + " at for child id " + childId + " (" + valueString + ") could not be parsed as a number."
                                );
                            }

                            covariatesX[childIndex][i] = newValue;
                            i++;

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

        // Adjust for covariates
        if (covariates.length > 0) {

            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();

            for (String phenoName : phenoNames) {

                double[] phenos = phenoMap.get(phenoName);

                int nValidValues = nValidValuesMap.get(phenoName);

                int[] index = new int[nValidValues];
                double[] y = new double[nValidValues];
                double[][] x = new double[nValidValues][covariatesX.length];

                int j = 0;
                for (int i = 0; i < phenos.length; i++) {

                    if (!Double.isNaN(phenos[i])) {

                        index[j] = i;
                        y[j] = phenos[i];
                        x[j] = covariatesX[i];
                        j++;

                    }
                }

                regression.newSampleData(y, x);

                try {

                    double[] residuals = regression.estimateResiduals();

                    double[] newPhenos = Arrays.copyOf(phenos, phenos.length);

                    for (int i = 0; i < index.length; i++) {

                        j = index[i];
                        newPhenos[j] = residuals[i];

                    }

                    phenoMap.put(phenoName, newPhenos);

                    for (int i = 0; i < phenos.length; i++) {

                        if (Double.isNaN(phenos[i]) && !Double.isNaN(newPhenos[i])) {

                            System.out.println(phenoName + ": " + phenos[i] + " - " + newPhenos[i]);
                            throw new IllegalArgumentException("NA issue");

                        }

                    }

                    System.out.println(phenoName + ": " + phenos[0] + "," + phenos[1] + "," + phenos[2] + "," + phenos[3] + "," + phenos[4] + "," + phenos[5]);
                    System.out.println(phenoName + ": " + newPhenos[0] + "," + newPhenos[1] + "," + newPhenos[2] + "," + newPhenos[3] + "," + newPhenos[4] + "," + newPhenos[5]);

                } catch (SingularMatrixException singularMatrixException) {

                    throw new IllegalArgumentException("Singular matrix obtained when adjusting " + phenoName + " for the given covariates.", singularMatrixException);

                }
            }
        }
    }
}
