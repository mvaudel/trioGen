package no.uib.triogen.model.simple_score;

import java.io.File;
import java.util.ArrayList;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class provides variants to include in scoring along with weights.
 *
 * @author Marc Vaudel
 */
public class VariantWeightList {

    /**
     * The variant list.
     */
    public final VariantList variantList;
    /**
     * The effect alleles.
     */
    public final String[] effectAllele;
    /**
     * The weights.
     */
    public final double[] weights;

    /**
     * Constructor.
     *
     * @param variantId the ids of the variant to include
     * @param chromosome the chromosome name where to look for
     * @param effectAllele the effect allele
     * @param position the position of the variant
     * @param weights the weights to use for the score
     */
    public VariantWeightList(
            String[] variantId,
            String[] chromosome,
            int[] position,
            String[] effectAllele,
            double[] weights
    ) {
        
        this.variantList = new VariantList(variantId, chromosome, position);
        this.effectAllele = effectAllele;
        this.weights = weights;
        
    }

    /**
     * Parses a variant weight list from a file. First lines starting with '#' are
     * ignored.
     *
     * @param variantFile the file
     *
     * @return an instance of variant weight list
     */
    public static VariantWeightList getVariantWeightList(
            File variantFile
    ) {

        ArrayList<String> variantIdList = new ArrayList<>();
        ArrayList<String> chromosomeList = new ArrayList<>();
        ArrayList<Integer> positionList = new ArrayList<>();
        ArrayList<String> effectAlleleList = new ArrayList<>();
        ArrayList<Double> weightList = new ArrayList<>();

        int lineNumber = 0;

        try ( SimpleFileReader reader = SimpleFileReader.getFileReader(variantFile)) {

            String line;
            while ((line = reader.readLine()) != null
                    && (line.trim().length() == 0 || line.charAt(0) == '#')) {

                lineNumber++;

            }

            while ((line = reader.readLine()) != null) {

                line = line.trim();

                if (line.length() > 0) {

                    lineNumber++;

                    String[] lineSplit = line.split(IoUtils.SEPARATOR);

                    if (lineSplit.length < 6) {

                        throw new IllegalArgumentException(
                                lineSplit.length + " elements found at line " + lineNumber + " where at least 6 expected. Please make sure that the file is tab-separated.\n" + line
                        );
                    }

                    String variantId = lineSplit[0];
                    String chromosome = lineSplit[1];
                    String positionString = lineSplit[2];
                    String effectAlleleString = lineSplit[3];
                    String weightString = lineSplit[4];

                    int position;
                    try {

                        position = Integer.parseInt(positionString);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("Position (" + positionString + ") could not be parsed as integer for variant " + variantId + " at line " + lineNumber + ".");

                    }

                    double weight;
                    try {

                        weight = Double.parseDouble(weightString);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("Weight (" + weightString + ") could not be parsed as number for variant " + variantId + " at line " + lineNumber + ".");

                    }

                    variantIdList.add(variantId);
                    chromosomeList.add(chromosome);
                    positionList.add(position);
                    effectAlleleList.add(effectAlleleString);
                    weightList.add(weight);

                }
            }
        }

        if (variantIdList.isEmpty()) {
            
            System.out.println("Warning: No variant found in " + variantFile + ". Please verify that the header is not commented by \"#\"");

        }

        return new VariantWeightList(
                variantIdList.toArray(new String[variantIdList.size()]),
                chromosomeList.toArray(new String[chromosomeList.size()]),
                positionList.stream()
                        .mapToInt(a -> a)
                        .toArray(),
                effectAlleleList.toArray(new String[effectAlleleList.size()]),
                weightList.stream()
                        .mapToDouble(a -> a)
                        .toArray()
        );
    }
}
