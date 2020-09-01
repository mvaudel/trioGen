package no.uib.triogen.model.covariates;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
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
    public final ConcurrentHashMap<String, int[]> indexMap;

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
        indexMap = new ConcurrentHashMap<>(phenotypeNames.length);

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

            int[] index = new int[nValidValues];
            double[] y = new double[nValidValues];
            double[][] x = new double[nValidValues][covariates.length + 1];

            // DEBUG
            boolean covariateNA = false;

            int j = 0;
            for (int i = 0; i < phenos.length; i++) {

                if (!Double.isNaN(phenos[i])) {

                    index[j] = i;
                    y[j] = phenos[i];

                    for (int k = 0; k < covariates.length; k++) {

                        String covariate = covariates[k];

                        double[] covariateValues = phenotypesHandler.phenoMap.get(covariate);

                        x[j][k] = covariateValues[i];

                        // DEBUG
                        if (Double.isNaN(covariateValues[i])) {
                            covariateNA = true;
                        }

                    }

                    x[j][covariates.length] = 1.0;

                    j++;

                }

            }

            // DEBUG
            if (covariateNA) {

                System.out.println("Covariate NA");

            }

            Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(x, false);

            SingularValueDecomposition svd = new SingularValueDecomposition(xMatrix);

            uMap.put(phenoName, svd.getU());
            utMap.put(phenoName, svd.getUT());
            rankMap.put(phenoName, svd.getRank());
            rawPhenoValues.put(phenoName, y);
            indexMap.put(phenoName, index);

        } else {

            double[] phenos = phenotypesHandler.phenoMap.get(phenoName);

            int nValidValues = phenotypesHandler.nValidValuesMap.get(phenoName);

            int[] index = new int[nValidValues];
            double[] y = new double[nValidValues];

            int j = 0;
            for (int i = 0; i < phenos.length; i++) {

                if (!Double.isNaN(phenos[i])) {

                    index[j] = i;
                    y[j] = phenos[i];

                    j++;

                }
            }

            rawPhenoValues.put(phenoName, y);

            indexMap.put(phenoName, index);

        }
    }

    /**
     * Removes the contribution of the covariates from the phenotypes and returns the new values.
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
