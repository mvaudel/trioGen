package no.uib.triogen.processing.mendelian_check;

import io.airlift.compress.zstd.ZstdDecompressor;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.genotypes.bgen.iterator.VariantIterator;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.model.mendelian_error.MendelianErrorEstimator;

/**
 * Runnable for the LD matrix writer.
 *
 * @author Marc Vaudel
 */
public class MendelianCheckRunnable implements Runnable {

    /**
     * The iterator.
     */
    private final VariantIterator iterator;
    /**
     * The index of the bgen file to process.
     */
    private final BgenIndex bgenIndex;
    /**
     * The reader for the bgen file to process.
     */
    private final BgenFileReader bgenFileReader;
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
     * The allele frequency threshold.
     */
    private final double alleleFrequencyThreshold;
    /**
     * The decompressor to use.
     */
    private final ZstdDecompressor decompressor = new ZstdDecompressor();

    /**
     * Constructor.
     *
     * @param writer The writer to use.
     * @param iterator The variant iterator.
     * @param bgenIndex The index of the bgen file.
     * @param bgenFileReader The reader for the bgen file.
     * @param childToParentMap The map of trios.
     * @param alleleFrequencyThreshold The allele frequency threshold.
     * values lower than threshold are not included (inclusive).
     * @param logger The logger.
     */
    public MendelianCheckRunnable(
            SimpleFileWriter writer,
            VariantIterator iterator,
            BgenIndex bgenIndex,
            BgenFileReader bgenFileReader,
            ChildToParentMap childToParentMap,
            double alleleFrequencyThreshold,
            SimpleCliLogger logger
    ) {

        this.writer = writer;
        this.iterator = iterator;
        this.bgenIndex = bgenIndex;
        this.bgenFileReader = bgenFileReader;
        this.childToParentMap = childToParentMap;
        this.alleleFrequencyThreshold = alleleFrequencyThreshold;
        this.logger = logger;

    }

    @Override
    public void run() {

        try {

            Integer tempIndex;
            while ((tempIndex = iterator.next()) != null && !canceled) {

                int variantIndex = tempIndex;
                VariantInformation variantInformation = bgenIndex.variantInformationArray[variantIndex];

                if (variantInformation.alleles.length > 1) {

                    BgenVariantData variantData = bgenFileReader.getVariantData(variantIndex);
                    variantData.parse(
                            childToParentMap,
                            decompressor
                    );

                    // Check if any allele passes the frequency threshold
                    int[] testedAlleleIndexes = IntStream.range(1, variantData.getOrderedAlleles().length)
                            .filter(
                                    alleleIndex -> variantData.getAlleleFrequency(alleleIndex) > alleleFrequencyThreshold 
                                            && variantData.getAlleleFrequency(alleleIndex) < 1.0 - alleleFrequencyThreshold
                            )
                            .toArray();

                    if (testedAlleleIndexes.length > 0) {

                        for (int alleleI : testedAlleleIndexes) {

                            int freq_hMnt_1 = 0;
                            int freq_hMnt_2 = 0;
                            int freq_hPnt_1 = 0;
                            int freq_hPnt_2 = 0;

                            for (String childId : childToParentMap.children) {

                                String motherId = childToParentMap.getMother(childId);
                                String fatherId = childToParentMap.getFather(childId);

                                double[] haplotypes = variantData.getHaplotypes(childId, motherId, fatherId, alleleI);

                                if (haplotypes[0] <= -0.5) {

                                    freq_hMnt_1 += 1;

                                }
                                if (haplotypes[0] >= 1.5) {

                                    freq_hMnt_2 += 1;

                                }
                                if (haplotypes[3] == -1) {

                                    freq_hPnt_1 += 1;

                                }
                                if (haplotypes[3] == 2) {

                                    freq_hPnt_2 += 1;

                                }
                            }

                            double alleleFrequency = variantData.getAlleleFrequency(alleleI);

                            double exp_hMnt_1 = (1 - alleleFrequency) * (1 - alleleFrequency) * alleleFrequency * childToParentMap.children.length; // 001*
                            double exp_hMnt_2 = alleleFrequency * alleleFrequency * (1 - alleleFrequency) * childToParentMap.children.length; // 110*
                            double exp_hPnt_1 = exp_hMnt_1 * childToParentMap.children.length; // *100
                            double exp_hPnt_2 = exp_hMnt_2 * childToParentMap.children.length; // *011

                            double p_hMnt_1 = ((double) freq_hMnt_1) / exp_hMnt_1;
                            double p_hMnt_2 = ((double) freq_hMnt_2) / exp_hMnt_2;
                            double p_hPnt_1 = ((double) freq_hPnt_1) / exp_hPnt_1;
                            double p_hPnt_2 = ((double) freq_hPnt_2) / exp_hPnt_2;

                            double prevalenceBefore = MendelianErrorEstimator.estimateMendelianErrorPrevalence(variantData, childToParentMap, alleleI);

                            variantData.swapChildrenAlleles();

                            double prevalenceAfter = MendelianErrorEstimator.estimateMendelianErrorPrevalence(variantData, childToParentMap, alleleI);

                            variantData.swapChildrenAlleles();

                            writer.writeLine(
                                    variantInformation.contig,
                                    Integer.toString(variantInformation.position),
                                    variantInformation.id,
                                    variantInformation.rsId,
                                    variantInformation.alleles[alleleI],
                                    variantInformation.getOtherAllele(alleleI),
                                    Double.toString(alleleFrequency),
                                    Double.toString(freq_hMnt_1),
                                    Double.toString(freq_hMnt_2),
                                    Double.toString(freq_hPnt_1),
                                    Double.toString(freq_hPnt_2),
                                    Double.toString(exp_hMnt_1),
                                    Double.toString(exp_hMnt_2),
                                    Double.toString(exp_hPnt_1),
                                    Double.toString(exp_hPnt_2),
                                    Double.toString(p_hMnt_1),
                                    Double.toString(p_hMnt_2),
                                    Double.toString(p_hPnt_1),
                                    Double.toString(p_hPnt_2),
                                    Double.toString(prevalenceBefore),
                                    Double.toString(prevalenceAfter)
                            );
                        }
                    }
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
}
