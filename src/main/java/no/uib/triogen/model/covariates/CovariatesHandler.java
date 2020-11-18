package no.uib.triogen.model.covariates;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.utils.Utils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Handler for covariates.
 *
 * @author Marc Vaudel
 */
public class CovariatesHandler {

    /**
     * The number of neighbors to use when masking missing covariates.
     */
    public final int nNeighbors = 200;
    /**
     * Phenotype name to covariates.
     */
    public final ConcurrentHashMap<String, String[]> covariatesMap;
    /**
     * The u for each phenotype.
     */
    public final ConcurrentHashMap<String, RealMatrix> uMap;
    /**
     * The transpose of u for each phenotype.
     */
    public final ConcurrentHashMap<String, RealMatrix> utMap;
    /**
     * The effective numerical rank for each phenotype.
     */
    public final ConcurrentHashMap<String, Integer> rankMap;
    /**
     * The phenotypes prior to adjustment for covariates.
     */
    public final ConcurrentHashMap<String, double[]> rawPhenoValues;
    /**
     * The index of the children in the raw pheno values.
     */
    public final ConcurrentHashMap<String, int[]> originalIndexMap;
    /**
     * The index of the children in the adjusted pheno values.
     */
    public final ConcurrentHashMap<String, int[]> adjustedIndexMap;

    /**
     * Constructor.
     *
     * @param phenotypeNames The names of the phenotypes to load.
     * @param phenotypesHandler The phenotype handler to use.
     * @param covariatesGeneral The array of covariates to use for all
     * phenotypes.
     * @param covariatesSpecific The map of covariates to use for phenotypes
     * specifically.
     */
    public CovariatesHandler(
            String[] phenotypeNames,
            PhenotypesHandler phenotypesHandler,
            String[] covariatesGeneral,
            HashMap<String, TreeSet<String>> covariatesSpecific
    ) {

        this.covariatesMap = new ConcurrentHashMap<>(phenotypeNames.length);
        uMap = new ConcurrentHashMap<>(phenotypeNames.length);
        utMap = new ConcurrentHashMap<>(phenotypeNames.length);
        rankMap = new ConcurrentHashMap<>(phenotypeNames.length);
        rawPhenoValues = new ConcurrentHashMap<>(phenotypeNames.length);
        originalIndexMap = new ConcurrentHashMap<>(phenotypeNames.length);
        adjustedIndexMap = new ConcurrentHashMap<>(phenotypeNames.length);

        Arrays.stream(phenotypeNames)
                .parallel()
                .forEach(
                        phenoName -> singularValueDecomposition(
                                phenoName,
                                phenotypesHandler,
                                covariatesGeneral,
                                covariatesSpecific
                        )
                );

    }

    /**
     * Runs svd for the given phenotype and fills the relevant maps.
     *
     * @param phenoName The name of the phenotype to process.
     * @param phenotypesHandler The phenotype handler to use.
     * @param covariatesGeneral The array of covariates to use for all
     * phenotypes.
     * @param covariatesSpecific The map of covariates to use for phenotypes
     * specifically.
     */
    private void singularValueDecomposition(
            String phenoName,
            PhenotypesHandler phenotypesHandler,
            String[] covariatesGeneral,
            HashMap<String, TreeSet<String>> covariatesSpecific
    ) {

        String[] covariates = Stream.concat(
                Arrays.stream(covariatesGeneral),
                covariatesSpecific.get(phenoName).stream()
        )
                .distinct()
                .toArray(String[]::new);

        covariatesMap.put(phenoName, covariates);

        if (covariates.length > 0) {

            double[] phenos = phenotypesHandler.phenoMap.get(phenoName);

            int nValidValues = phenotypesHandler.nValidValuesMap.get(phenoName);

            int[] originalIndex = new int[nValidValues];
            int[] adjustedIndex = new int[phenos.length];
            Arrays.fill(adjustedIndex, -1);
            
            double[] y = new double[nValidValues];
            double[][] x = new double[nValidValues][covariates.length + 1];

            HashMap<Integer, ArrayList<Integer>> naKtoJ = new HashMap<>();

            int j = 0;
            for (int i = 0; i < phenos.length; i++) {

                if (!Double.isNaN(phenos[i])) {

                    originalIndex[j] = i;
                    adjustedIndex[i] = j;
                    
                    y[j] = phenos[i];

                    for (int k = 0; k < covariates.length; k++) {

                        String covariate = covariates[k];

                        double[] covariateValues = phenotypesHandler.phenoMap.get(covariate);

                        double covariateValue = covariateValues[i];

                        if (!Double.isNaN(covariateValue)) {

                            x[j][k] = covariateValue;

                        } else {

                            ArrayList<Integer> naJ = naKtoJ.get(k);

                            if (naJ == null) {

                                naJ = new ArrayList<>();
                                naKtoJ.put(k, naJ);

                            }

                            naJ.add(j);

                        }
                    }

                    x[j][covariates.length] = 1.0;

                    j++;

                }
            }

            // Mask missing covariates
            if (!naKtoJ.isEmpty()) {

                for (Entry<Integer, ArrayList<Integer>> entry : naKtoJ.entrySet()) {

                    int k = entry.getKey();

                    TreeMap<Double, ArrayList<Double>> phenoToCovariateValueMap = new TreeMap<>();

                    String covariate = covariates[k];

                    double[] covariateValues = phenotypesHandler.phenoMap.get(covariate);

                    for (int i = 0; i < phenos.length; i++) {

                        double phenoValue = phenos[i];
                        double covariateValue = covariateValues[i];

                        if (!Double.isNaN(phenoValue) && !Double.isNaN(covariateValue)) {

                            ArrayList<Double> covariatesAtValue = phenoToCovariateValueMap.get(phenoValue);

                            if (covariatesAtValue == null) {

                                covariatesAtValue = new ArrayList<>(1);
                                phenoToCovariateValueMap.put(phenoValue, covariatesAtValue);

                            }

                            covariatesAtValue.add(covariateValue);

                        }
                    }

                    HashMap<Double, Integer> covariatePhenoIndex = new HashMap<>(y.length);

                    for (int i = 0; i < y.length; i++) {

                        covariatePhenoIndex.put(y[i], i);

                    }

                    for (int i : entry.getValue()) {

                        ArrayList<Double> neighbors = new ArrayList<>(nNeighbors);

                        double phenoValue = y[i];
                        int phenoIndex = covariatePhenoIndex.get(phenoValue);

                        ArrayList<Double> covariatesAtPheno = phenoToCovariateValueMap.get(phenoValue);
                        int covariatesAtPhenoSize = 0;

                        if (covariatesAtPheno != null) {

                            neighbors.addAll(covariatesAtPheno);
                            covariatesAtPhenoSize = covariatesAtPheno.size();

                        }

                        int toFill = (nNeighbors - covariatesAtPhenoSize) / 2;
                        int filled = 0;

                        for (j = phenoIndex - 1; j >= 0 && filled < toFill; j--) {

                            covariatesAtPheno = phenoToCovariateValueMap.get(y[j]);

                            if (covariatesAtPheno != null) {

                                filled += covariatesAtPheno.size();

                                neighbors.addAll(covariatesAtPheno);

                            }
                        }

                        filled = 0;

                        for (j = phenoIndex + 1; j < y.length && filled < toFill; j--) {

                            covariatesAtPheno = phenoToCovariateValueMap.get(y[j]);

                            if (covariatesAtPheno != null) {

                                filled += covariatesAtPheno.size();

                                neighbors.addAll(covariatesAtPheno);

                            }
                        }

                        Collections.sort(neighbors);

                        double maskValue = Utils.medianSorted(neighbors);

                        x[i][k] = maskValue;

                    }
                }
            }

            Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(x, false);

            SingularValueDecomposition svd = new SingularValueDecomposition(xMatrix);

            uMap.put(phenoName, svd.getU());
            utMap.put(phenoName, svd.getUT());
            rankMap.put(phenoName, svd.getRank());
            rawPhenoValues.put(phenoName, y);
            originalIndexMap.put(phenoName, originalIndex);
            adjustedIndexMap.put(phenoName, adjustedIndex);

        } else {

            double[] phenos = phenotypesHandler.phenoMap.get(phenoName);

            int nValidValues = phenotypesHandler.nValidValuesMap.get(phenoName);

            int[] originalIndex = new int[nValidValues];
            int[] adjustedIndex = new int[phenos.length];
            Arrays.fill(adjustedIndex, -1);
            
            double[] y = new double[nValidValues];

            int j = 0;
            for (int i = 0; i < phenos.length; i++) {

                if (!Double.isNaN(phenos[i])) {

                    originalIndex[j] = i;
                    adjustedIndex[i] = j;
                    y[j] = phenos[i];

                    j++;

                }
            }

            rawPhenoValues.put(phenoName, y);

            originalIndexMap.put(phenoName, originalIndex);
            adjustedIndexMap.put(phenoName, adjustedIndex);

        }
    }

    /**
     * Removes the contribution of the covariates from the phenotypes and
     * returns the new values.
     *
     * @param phenoName The name of the phenotype.
     *
     * @return The values adjusted for covariates.
     */
    public double[] getAdjustedValues(
            String phenoName
    ) {

        double[] values = rawPhenoValues.get(phenoName);

        RealMatrix utMatrix = utMap.get(phenoName);

        if (utMatrix != null) {

            double[] covariateBaseValues = utMatrix.operate(values);

            for (int i = rankMap.get(phenoName) + 1; i < covariateBaseValues.length; i++) {

                covariateBaseValues[i] = 0.0;

            }

            RealMatrix uMatrix = uMap.get(phenoName);
            double[] covariatesContribution = uMatrix.operate(covariateBaseValues);

            for (int i = 0; i < values.length; i++) {

                values[i] = values[i] - covariatesContribution[i];

            }
        }

        return values;

    }

    /**
     * Subtracts the contribution of the covariates to the columns in x.
     *
     * @param phenoName The name of the phenotype.
     * @param x The values to project.
     *
     * @return The projected values.
     */
    public double[][] getAdjustedValues(
            String phenoName,
            double[][] x
    ) {

        RealMatrix utMatrix = utMap.get(phenoName);

        if (utMatrix != null) {

            Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(x, false);

            RealMatrix covariateBaseValues = utMatrix.multiply(xMatrix);

            for (int i = rankMap.get(phenoName) + 1; i < covariateBaseValues.getRowDimension(); i++) {

                for (int j = 0; j < covariateBaseValues.getColumnDimension(); j++) {

                    covariateBaseValues.setEntry(i, j, 0.0);

                }
            }

            RealMatrix uMatrix = uMap.get(phenoName);
            RealMatrix covariatesContributionMatrix = uMatrix.multiply(covariateBaseValues);

            RealMatrix resultMatrix = xMatrix.subtract(covariatesContributionMatrix);

            return resultMatrix.getData();

        } else {

            return x;

        }
    }

    /**
     * Removes the raw pheno values to save memory.
     */
    public void trim() {

        rawPhenoValues.clear();

    }
}
