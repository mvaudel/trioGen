package no.uib.triogen.processing.prs;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.model.trio_genotypes.VariantList;
import static no.uib.triogen.processing.prs.PrsPruner.VARIABLE_WILDCARD;

/**
 * Utils for PRS calculation.
 *
 * @author Marc Vaudel
 */
public class PrsUtils {

    /**
     * The modes of scoring
     */
    public enum ScoringMode {

        lead(0),
        weighted(1);

        private ScoringMode(int index) {

            this.index = index;

        }

        public final int index;

    }

    /**
     * For each variant, parses the weight and effect size.
     *
     * @param trainingFile
     * @param betaColumnPattern
     * @param seColumnPattern
     * @param variableNames The names of the variables
     * @param variantList
     * @param pValueThreshold
     * @param scoringMode
     *
     * @return The scoring data in a map, chromosome to lead variant to scoring
     * variant to effect allele to weight and effect size.
     */
    public static HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> parseScoringData(
            File trainingFile,
            String betaColumnPattern,
            String seColumnPattern,
            String[] variableNames,
            VariantList variantList,
            double pValueThreshold,
            ScoringMode scoringMode
    ) {

        HashMap<String, HashMap<String, HashMap<String, HashMap<String, double[]>>>> scoringData = new HashMap<>();

        int[] betaColumnIndexes = new int[variableNames.length];
        Arrays.fill(betaColumnIndexes, -1);

        int[] seColumnIndexes = new int[variableNames.length];
        Arrays.fill(seColumnIndexes, -1);

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(trainingFile)) {

            String line = reader.readLine();

            String[] lineSplit = line.split(SEPARATOR);

            for (int i = 0; i < lineSplit.length; i++) {

                for (int j = 0; j < variableNames.length; j++) {

                    String variable = variableNames[j];

                    String betaColumn = betaColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(betaColumn)) {

                        betaColumnIndexes[j] = i;

                    }

                    String seColumn = seColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(seColumn)) {

                        seColumnIndexes[j] = i;

                    }
                }
            }

            for (int j = 0; j < variableNames.length; j++) {

                String variable = variableNames[j];

                if (betaColumnIndexes[j] == -1) {

                    String betaColumn = betaColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);
                    throw new IllegalArgumentException("Effect size column '" + betaColumn + "' not found for '" + variable + "'.");

                }

                if (seColumnIndexes[j] == -1) {

                    String seColumn = seColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);
                    throw new IllegalArgumentException("Standard error column '" + seColumn + "' not found for '" + variable + "'.");

                }
            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(SEPARATOR);

                String leadVariantId = lineSplit[0];
                double leadP = Double.parseDouble(lineSplit[2]);

                String variantId = lineSplit[3];
                String variantRsid = lineSplit[4];

                if (variantList != null && variantList.contains(variantId)
                        || variantList != null && variantList.contains(variantRsid)
                        || variantList == null && leadP <= pValueThreshold) {

                    if (scoringMode == ScoringMode.weighted || leadVariantId.equals(variantId)) {

                        String chromosome = lineSplit[5];
                        String ea = lineSplit[8];

                        double[] result = new double[2 * variableNames.length];

                        for (int j = 0; j < variableNames.length; j++) {

                            String betaString = lineSplit[betaColumnIndexes[j]];

                            double beta = Double.parseDouble(betaString);

                            result[j] = beta;

                            String seString = lineSplit[seColumnIndexes[j]];

                            double se = Double.parseDouble(seString);

                            result[j + variableNames.length] = se;

                        }

                        HashMap<String, HashMap<String, HashMap<String, double[]>>> chromosomeValues = scoringData.get(chromosome);

                        if (chromosomeValues == null) {

                            chromosomeValues = new HashMap<>();
                            scoringData.put(chromosome, chromosomeValues);

                        }

                        HashMap<String, HashMap<String, double[]>> leadVariantValues = chromosomeValues.get(leadVariantId);

                        if (leadVariantValues == null) {

                            leadVariantValues = new HashMap<>(1);
                            chromosomeValues.put(leadVariantId, leadVariantValues);

                        }

                        HashMap<String, double[]> variantValues = leadVariantValues.get(variantId);

                        if (variantValues == null) {

                            variantValues = new HashMap<>(1);
                            leadVariantValues.put(variantId, variantValues);

                        }

                        variantValues.put(ea, result);

                    }
                }
            }
        }

        return scoringData;

    }
}
