package no.uib.triogen.io.genotypes.iterators;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.util.ArrayList;
import no.uib.triogen.io.genotypes.GenotypesIterator;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import static no.uib.triogen.io.genotypes.vcf.generic.VcfIterator.getVcfIndexFile;
import no.uib.triogen.io.genotypes.vcf.generic.VcfVariant;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.maf.MafEstimator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Iterator for a single SNP and the SNPs around in a given BP window.
 *
 * @author Marc Vaudel
 */
public class TargetGenotypesIterator implements GenotypesIterator {

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
     * Window start for the iterator.
     */
    public final int windowStart;
    /**
     * Window end for the iterator.
     */
    public final int windowEnd;
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
     * Genotypes in the window.
     */
    private ArrayList<GenotypesProvider> window = new ArrayList<>();
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

        windowStart = targetSnpBpStart - downStreamDistance;
        windowEnd = targetSnpBpEnd + upStreamDistance;

        this.mafThreshold = mafThreshold;
        this.childToParentMap = childToParentMap;

    }

    private void init() {

        File indexFile = getVcfIndexFile(vcfFile);
        VCFFileReader reader = new VCFFileReader(vcfFile, indexFile);
        CloseableIterator<VariantContext> iterator = reader.query(targetSnpContig,
                windowStart,
                windowEnd
        );

        VariantContext variantContext;

        while ((variantContext = iterator.next()) != null) {

            GenotypesProvider tempGenotypesProvider = new VcfVariant(
                    variantContext,
                    false
            );

            tempGenotypesProvider.parse();

            double maf = MafEstimator.getMaf(tempGenotypesProvider,
                    childToParentMap
            );

            if (maf >= mafThreshold) {

                tempGenotypesProvider.setParentP0s(childToParentMap.children, childToParentMap);

                window.add(tempGenotypesProvider);

                String vcfVariantId = variantContext.getID();

                if (vcfVariantId.equals(targetSnpId)) {

                    genotypesProvider = new VcfVariant(
                            variantContext,
                            true
                    );
                }
            }
        }

        iterator.close();
        reader.close();

    }

    @Override
    public GenotypesProvider next() {

        init();

        semaphore.acquire();

        GenotypesProvider result = genotypesProvider;

        genotypesProvider = null;

        semaphore.release();

        return result;

    }

    @Override
    public GenotypesProvider[] getGenotypesInRange(
            String contig,
            int startBp,
            int endBp
    ) {

        if (startBp < windowStart) {

            throw new IndexOutOfBoundsException("window start (" + startBp + ") out of range (" + startBp + " - " + endBp + ").");

        }

        if (endBp < windowEnd) {

            throw new IndexOutOfBoundsException("window start (" + startBp + ") out of range (" + startBp + " - " + endBp + ").");

        }

        return window.stream()
                .filter(
                        result -> result.getBp() >= startBp && result.getBp() <= endBp
                )
                .toArray(
                        GenotypesProvider[]::new
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
