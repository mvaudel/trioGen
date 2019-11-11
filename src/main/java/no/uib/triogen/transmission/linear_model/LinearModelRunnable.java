package no.uib.triogen.transmission.linear_model;

import java.util.Arrays;
import java.util.HashMap;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.TransmissionModel;
import no.uib.triogen.model.pheno.PhenotypesHandler;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 * Runnable for the linear model association.
 *
 * @author Marc Vaudel
 */
public class LinearModelRunnable implements Runnable {

    /**
     * The variants iterator.
     */
    private final VariantIterator iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The handler for phenotypes.
     */
    private final PhenotypesHandler phenotypesHandler;
    /**
     * The output writer.
     */
    private final SimpleFileWriter outputWriter;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;

    /**
     * Constructor.
     *
     * @param iterator the variants iterator
     * @param childToParentMap the child to parent map
     * @param phenotypesHandler the phenotypes handler
     * @param outputWriter the output writer
     */
    public LinearModelRunnable(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            PhenotypesHandler phenotypesHandler,
            SimpleFileWriter outputWriter
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.phenotypesHandler = phenotypesHandler;
        this.outputWriter = outputWriter;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider tempGenotypesProvider;
            while ((tempGenotypesProvider = iterator.next()) != null && !canceled) {

                GenotypesProvider genotypesProvider = tempGenotypesProvider;
                Arrays.stream(TransmissionModel.values())
                        .parallel()
                        .forEach(
                                transmissionModel -> runLinearModel(
                                        transmissionModel,
                                        genotypesProvider
                                )
                        );
            }

        } catch (Throwable t) {

            canceled = true;
            t.printStackTrace();

        }
    }

    /**
     * Runs the linear model for the given transmission model.
     * 
     * @param transmissionModel the model to use for transmission
     * @param genotypesProvider the genotypes provider
     */
    public void runLinearModel(
            TransmissionModel transmissionModel,
            GenotypesProvider genotypesProvider
    ) {

        phenotypesHandler.phenoMap.entrySet()
                .parallelStream()
                .forEach(
                        entry -> runLinearModel(
                                transmissionModel,
                                genotypesProvider,
                                entry.getKey(),
                                entry.getValue()
                        )
                );
    }

    /**
     * Runs the linear model for a transmission model and a phenotype name.
     * 
     * @param transmissionModel the model to use for transmission
     * @param genotypesProvider the genotypes provider
     * @param phenoName the phenotype name
     * @param phenotypes the phenotype values
     */
    public void runLinearModel(
            TransmissionModel transmissionModel,
            GenotypesProvider genotypesProvider,
            String phenoName,
            HashMap<String, Double> phenotypes
    ) {

        SimpleRegression simpleRegression = new SimpleRegression();

        for (String childId : childToParentMap.children) {

            double[] hs = genotypesProvider.getH(childToParentMap, childId);
            double x = transmissionModel.getRegressionValue(hs);

            if (!Double.isNaN(x) && !Double.isInfinite(x)) {

                double y = phenotypes.get(childId);

                if (!Double.isNaN(y) && !Double.isInfinite(y)) {

                    simpleRegression.addData(x, y);

                }
            }
        }

        String line = String.join(
                Utils.separator,
                phenoName,
                genotypesProvider.getVariantID(),
                transmissionModel.name(),
                transmissionModel.h,
                Double.toString(simpleRegression.getSlope()),
                Double.toString(simpleRegression.getSlopeStdErr()),
                Double.toString(simpleRegression.getSlopeConfidenceInterval()),
                Double.toString(simpleRegression.getRSquare()),
                Double.toString(simpleRegression.getSignificance()),
                Long.toString(simpleRegression.getN())
        );
        outputWriter.writeLine(line);

    }
}
