package no.uib.triogen.transmission.linear_model;

import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.TransmissionModel;
import no.uib.triogen.model.pheno.PhenotypesHandler;

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

            GenotypesProvider genotypesProvider;
            while ((genotypesProvider = iterator.next()) != null && !canceled) {

                

            }
            
        } catch (Throwable t) {

            canceled = true;
            t.printStackTrace();

        }
    }
    
    public void runLinearModel(
            TransmissionModel transmissionModel
    ) {
        
        
    }

}
