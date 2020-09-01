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
            double[][] x = new double[nValidValues][nValidValues];

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

            utMap.put(phenoName, svd.getUT());

            rankMap.put(phenoName, svd.getRank());

            rawPhenoValues.put(phenoName, y);

            indexMap.put(phenoName, index);

            // DEBUG
            File degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_u.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                for (int i = 0; i < svd.getU().getRowDimension(); i++) {

                    String line = Arrays.stream(svd.getU().getRow(i))
                            .mapToObj(
                                    value -> Double.toString(value)
                            )
                            .collect(
                                    Collectors.joining("\t")
                            );
                    writer.writeLine(line);

                }
            }

            degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_s.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                for (int i = 0; i < svd.getS().getRowDimension(); i++) {

                    String line = Arrays.stream(svd.getS().getRow(i))
                            .mapToObj(
                                    value -> Double.toString(value)
                            )
                            .collect(
                                    Collectors.joining("\t")
                            );
                    writer.writeLine(line);

                }
            }

            degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_v.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                for (int i = 0; i < svd.getV().getRowDimension(); i++) {

                    String line = Arrays.stream(svd.getV().getRow(i))
                            .mapToObj(
                                    value -> Double.toString(value)
                            )
                            .collect(
                                    Collectors.joining("\t")
                            );
                    writer.writeLine(line);

                }
            }

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
     * Projects the phenotypes of the given name orthogonally to the covariates,
     * ie returns ut*values on the dimensions higher than the numerical rank.
     *
     * @param phenoName The name of the phenotype.
     *
     * @return The projected values.
     */
    public double[] getProjectedValues(
            String phenoName
    ) {

        double[] values = rawPhenoValues.get(phenoName);

        RealMatrix utMatrix = utMap.get(phenoName);

        if (utMatrix != null) {

            double[] product = utMatrix.operate(values);

            int resultLength = product.length - rankMap.get(phenoName);

            double[] result = new double[resultLength];

            System.arraycopy(product, resultLength, result, 0, resultLength);

            // DEBUG
            File degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_rawPheno.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                String line = Arrays.stream(values)
                        .mapToObj(
                                value -> Double.toString(value)
                        )
                        .collect(
                                Collectors.joining("\t")
                        );
                writer.writeLine(line);
            }

            degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_product.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                String line = Arrays.stream(product)
                        .mapToObj(
                                value -> Double.toString(value)
                        )
                        .collect(
                                Collectors.joining("\t")
                        );
                writer.writeLine(line);
            }

            degugFile = new File("/mnt/cargo/marc/triogen/svd/" + phenoName + "_result.gz");
            try ( SimpleFileWriter writer = new SimpleFileWriter(degugFile, true)) {

                String line = Arrays.stream(result)
                        .mapToObj(
                                value -> Double.toString(value)
                        )
                        .collect(
                                Collectors.joining("\t")
                        );
                writer.writeLine(line);
            }

            return result;

        } else {

            return values;

        }
    }

    /**
     * Projects the columns of the given values orthogonally to the covariates,
     * ie returns ut*values on the dimensions higher than the numerical rank.
     *
     * @param phenoName The name of the phenotype.
     * @param x The values to project.
     *
     * @return The projected values.
     */
    public double[][] getProjectedValues(
            String phenoName,
            double[][] x
    ) {

        RealMatrix utMatrix = utMap.get(phenoName);

        if (utMatrix != null) {

            Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(x, false);

            RealMatrix productMatrix = utMatrix.multiply(xMatrix);

            RealMatrix projectionMatrix = productMatrix.getSubMatrix(
                    rankMap.get(phenoName),
                    productMatrix.getRowDimension(),
                    0,
                    productMatrix.getColumnDimension()
            );

            return projectionMatrix.getData();

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