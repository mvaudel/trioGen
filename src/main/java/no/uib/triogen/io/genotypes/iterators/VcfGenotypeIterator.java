package no.uib.triogen.io.genotypes.iterators;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import static no.uib.triogen.io.genotypes.vcf.generic.VcfIterator.getVcfIndexFile;
import no.uib.triogen.io.genotypes.vcf.generic.VcfVariant;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.ChromosomeUtils;
import no.uib.triogen.model.maf.MafEstimator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Iterator based on a vcf file.
 *
 * @author Marc Vaudel
 */
public class VcfGenotypeIterator implements VariantIterator {

    /**
     * The iterator for the vcf file.
     */
    private final CloseableIterator<VariantContext> iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included (inclusive).
     */
    private final double mafThreshold;
    /**
     * The number of variants iterated.
     */
    private int nVariants = 0;
    /**
     * The coordinate of the start of the iteration.
     */
    private final int bpStart;
    /**
     * The coordinate of the end of the iteration.
     */
    private final int bpEnd;
    /**
     * Boolean indicating whether the iteration is finished.
     */
    private boolean finished = false;
    /**
     * Semaphore for the next method.
     */
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param vcfFile The vcf file to iterate.
     * @param contig The contig to iterate.
     * @param bpStart The start bp.
     * @param bpEnd The end bp.
     * @param mafThreshold The maf threshold to use.
     * @param childToParentMap The child to parent map to use.
     */
    public VcfGenotypeIterator(
            File vcfFile,
            String contig,
            int bpStart,
            int bpEnd,
            double mafThreshold,
            ChildToParentMap childToParentMap
    ) {

        this.mafThreshold = mafThreshold;
        this.childToParentMap = childToParentMap;
        this.bpStart = bpStart;
        this.bpEnd = bpEnd;

        int windowStart = Math.max(bpStart, 1);
        
        Integer chrLength = ChromosomeUtils.chromosomeLength37.get(contig);
        
        int windowEnd = chrLength != null && chrLength < bpEnd ? chrLength : bpEnd;

        File indexFile = getVcfIndexFile(vcfFile);
        VCFFileReader reader = new VCFFileReader(vcfFile, indexFile);
        iterator = reader.query(
                contig,
                windowStart,
                windowEnd
        );

    }

    @Override
    public GenotypesProvider next() {

        semaphore.acquire();

        VariantContext variantContext;

        while ((variantContext = iterator.next()) != null) {

            if (variantContext.getEnd() > bpEnd) {

                break;

            }

            if (variantContext.getStart() >= bpStart) {

                GenotypesProvider genotypesProvider = new VcfVariant(
                        variantContext
                );

                genotypesProvider.parse(childToParentMap);

                double maf = MafEstimator.getMaf(
                        genotypesProvider,
                        childToParentMap
                );

                if (maf >= mafThreshold) {

                    genotypesProvider.setParentP0s(childToParentMap.children, childToParentMap);

                    nVariants++;

                    semaphore.release();

                    return genotypesProvider;

                }
            }
        }

        finished = true;

        semaphore.release();

        return null;

    }

    @Override
    public int getnVariants() {

        return nVariants;

    }

    @Override
    public void close() {

        iterator.close();

    }

    @Override
    public boolean isFinished() {

        return finished;

    }

}
