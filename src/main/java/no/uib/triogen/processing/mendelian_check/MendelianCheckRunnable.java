package no.uib.triogen.processing.mendelian_check;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import static no.uib.triogen.cmd.ld.LdMatrixOptions.hardCalls;
import static no.uib.triogen.cmd.ld.LdMatrixOptions.maxDistance;
import static no.uib.triogen.cmd.ld.LdMatrixOptions.minR2;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.maf.MafEstimator;
import no.uib.triogen.model.mendelian_error.MendelianErrorEstimator;
import no.uib.triogen.model.trio_genotypes.VariantIndex;

/**
 * Runnable for the LD matrix writer.
 *
 * @author Marc Vaudel
 */
public class MendelianCheckRunnable implements Runnable, AutoCloseable {

    /**
     * The iterator.
     */
    private final VariantIterator iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;
    /**
     * Writer for the output.
     */
    private final SimpleFileWriter writer;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included (inclusive).
     */
    private final double mafThreshold;

    /**
     * Constructor.
     *
     * @param writer The writer to use.
     * @param iterator The variant iterator.
     * @param childToParentMap The map of trios.
     * @param mafThreshold The maf threshold. maf is computed in parents and
     * values lower than threshold are not included (inclusive).
     * @param logger The logger.
     */
    public MendelianCheckRunnable(
            SimpleFileWriter writer,
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            double mafThreshold,
            SimpleCliLogger logger
    ) {

        this.writer = writer;
        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.mafThreshold = mafThreshold;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider genotypesProvider;
            while ((genotypesProvider = iterator.next()) != null && !canceled) {
                
                genotypesProvider.parse(childToParentMap);

                double maf = MafEstimator.getMaf(genotypesProvider, childToParentMap);

                if (maf >= mafThreshold) {

                    double freq_h2_1 = 0.0;
                    double freq_h2_2 = 0.0;
                    double freq_h4_1 = 0.0;
                    double freq_h4_2 = 0.0;

                    for (String childId : childToParentMap.children) {

                        String motherId = childToParentMap.getMother(childId);
                        String fatherId = childToParentMap.getFather(childId);

                        short[] hs = genotypesProvider.getNAltH(childId, motherId, fatherId);

                        if (hs[1] == -1) {

                            freq_h2_1 += 1.0;

                        } 
                        if (hs[1] == 2) {

                            freq_h2_2 += 1.0;

                        } 
                        if (hs[3] == -1) {

                            freq_h4_1 += 1.0;

                        } 
                        if (hs[3] == 2) {

                            freq_h4_2 += 1.0;

                        }
                    }

                    freq_h2_1 /= childToParentMap.children.length;
                    freq_h2_2 /= childToParentMap.children.length;
                    freq_h4_1 /= childToParentMap.children.length;
                    freq_h4_2 /= childToParentMap.children.length;

                    double exp_h2_1 = (1 - maf) * (1 - maf) * maf; // 001*
                    double exp_h2_2 = maf * maf * (1 - maf); // 110*
                    double exp_h4_1 = exp_h2_1; // *100
                    double exp_h4_2 = exp_h2_2; // *011

                    double p_h2_1 = freq_h2_1 / exp_h2_1;
                    double p_h2_2 = freq_h2_2 / exp_h2_2;
                    double p_h4_1 = freq_h4_1 / exp_h4_1;
                    double p_h4_2 = freq_h4_2 / exp_h4_2;
                    
                    double prevalence = MendelianErrorEstimator.estimateMendelianErrorPrevalence(genotypesProvider, childToParentMap);

                    writer.writeLine(
                            genotypesProvider.getContig(),
                            Integer.toString(genotypesProvider.getBp()),
                            genotypesProvider.getVariantID(),
                            genotypesProvider.getRef(),
                            genotypesProvider.getAlt(),
                            genotypesProvider.genotyped() ? "1" : "0",
                            Double.toString(maf),
                            Double.toString(freq_h2_1),
                            Double.toString(freq_h2_2),
                            Double.toString(freq_h4_1),
                            Double.toString(freq_h4_2),
                            Double.toString(exp_h2_1),
                            Double.toString(exp_h2_2),
                            Double.toString(exp_h4_1),
                            Double.toString(exp_h4_2),
                            Double.toString(p_h2_1),
                            Double.toString(p_h2_2),
                            Double.toString(p_h4_1),
                            Double.toString(p_h4_2),
                            Double.toString(prevalence)
                    );
                }
            }

        } catch (Throwable t) {

            canceled = true;

            logger.logError(
                    Arrays.stream(t.getStackTrace())
                            .map(
                                    element -> element.toString()
                            )
                            .collect(Collectors.joining(" "))
            );

            t.printStackTrace();

        }
    }

    @Override
    public void close() throws Exception {

    }
}
