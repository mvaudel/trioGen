package no.uib.triogen.io.genotypes.iterators.window;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import static no.uib.triogen.io.genotypes.vcf.generic.VcfIterator.getVcfIndexFile;
import no.uib.triogen.io.genotypes.vcf.generic.VcfVariant;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.utils.SimpleSemaphore;
import no.uib.triogen.io.genotypes.WindowGenotypesIterator;
import no.uib.triogen.io.genotypes.iterators.VcfGenotypeIterator;
import no.uib.triogen.model.genome.ChromosomeUtils;

/**
 * Iterator for a single SNP and the SNPs around in a given BP window.
 *
 * @author Marc Vaudel
 */
public class TargetGenotypesIterator implements WindowGenotypesIterator {

    /**
     * The vcf file.
     */
    public final File vcfFile;
    /**
     * The contig of the target snp.
     */
    public final String targetSnpContig;
    /**
     * The id of the target snp.
     */
    public final String targetSnpId;
    /**
     * The bp start of the target snp.
     */
    public final int targetSnpBpStart;
    /**
     * The bp end of the target snp.
     */
    public final int targetSnpBpEnd;
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
     * Genotypes for the target snp.
     */
    private GenotypesProvider genotypesProvider = null;
    /**
     * Semaphore to access the target snp.
     */
    private SimpleSemaphore semaphore = new SimpleSemaphore(1);

    /**
     * Constructor.
     *
     * @param vcfFile The vcf file to get the genotypes from.
     * @param targetSnpId The id of target snp.
     * @param targetSnpContig The contig of the target snp.
     * @param targetSnpBpStart The start bp of the target snp.
     * @param targetSnpBpEnd The end bp of the target snp.
     * @param upStreamDistance The distance to load upstream the target snp.
     * @param downStreamDistance The distance to load downstream the target snp.
     * @param mafThreshold The maf threshold. maf is computed in parents and
     * values lower than threshold are not included (inclusive).
     * @param childToParentMap The map of trios.
     */
    public TargetGenotypesIterator(
            File vcfFile,
            String targetSnpId,
            String targetSnpContig,
            int targetSnpBpStart,
            int targetSnpBpEnd,
            int upStreamDistance,
            int downStreamDistance,
            double mafThreshold,
            ChildToParentMap childToParentMap
    ) {

        this.vcfFile = vcfFile;
        this.targetSnpContig = targetSnpContig;
        this.targetSnpId = targetSnpId;
        this.targetSnpBpStart = targetSnpBpStart;
        this.targetSnpBpEnd = targetSnpBpEnd;
        this.mafThreshold = mafThreshold;
        this.childToParentMap = childToParentMap;

        init();

    }

    /**
     * Sets up the iterator.
     */
    private void init() {

        int windowStart = Math.max(targetSnpBpStart, 1);
        
        Integer chrLength = ChromosomeUtils.chromosomeLength37.get(targetSnpContig);
        
        int windowEnd = chrLength != null && chrLength < targetSnpBpEnd ? chrLength : targetSnpBpEnd;

        File indexFile = getVcfIndexFile(vcfFile);
        VCFFileReader reader = new VCFFileReader(vcfFile, indexFile);
        CloseableIterator<VariantContext> iterator = reader.query(targetSnpContig,
                windowStart,
                windowEnd
        );

        VariantContext variantContext;

        while ((variantContext = iterator.next()) != null) {

            String vcfVariantId = variantContext.getID();

            if (vcfVariantId.equals(targetSnpId)) {

                genotypesProvider = new VcfVariant(variantContext, true);
                genotypesProvider.parse();
                genotypesProvider.setParentP0s(childToParentMap.children, childToParentMap);

            }
        }

        iterator.close();

        reader.close();

    }

    @Override
    public GenotypesProvider next() {

        semaphore.acquire();

        GenotypesProvider result = genotypesProvider;

        genotypesProvider = null;

        semaphore.release();

        return result;

    }

    @Override
    public VariantIterator getGenotypesInRange(
            String contig,
            int startBp,
            int endBp
    ) {

        return new VcfGenotypeIterator(
                vcfFile,
                contig,
                startBp,
                endBp,
                mafThreshold,
                childToParentMap
        );
    }

    @Override
    public void releaseMinBp(
            String contig,
            int minBp
    ) {

        // Nothing to do for this iterator.
    }

}
