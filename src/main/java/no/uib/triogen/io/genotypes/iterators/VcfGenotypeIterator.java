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
    private int nVariants = 0;
    private boolean finished = false;
    private final SimpleSemaphore semaphore = new SimpleSemaphore(1);

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

        File indexFile = getVcfIndexFile(vcfFile);
        VCFFileReader reader = new VCFFileReader(vcfFile, indexFile);
        iterator = reader.query(
                contig,
                bpStart,
                bpEnd
        );

    }

    @Override
    public GenotypesProvider next() {

        semaphore.acquire();

        VariantContext variantContext;

        while ((variantContext = iterator.next()) != null) {

            GenotypesProvider genotypesProvider = new VcfVariant(
                    variantContext,
                    false
            );

            genotypesProvider.parse();

            double maf = MafEstimator.getMaf(genotypesProvider,
                    childToParentMap
            );

            if (maf >= mafThreshold) {

                genotypesProvider.setParentP0s(childToParentMap.children, childToParentMap);

                nVariants++;

                semaphore.release();

                return genotypesProvider;

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
