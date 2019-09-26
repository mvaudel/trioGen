package no.uib.triogen.transmission.extraction;

import java.util.stream.Collectors;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.vcf.VcfIterator;
import no.uib.triogen.io.vcf.VcfLine;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Runnable for the extraction of transmitted alleles.
 *
 * @author Marc Vaudel
 */
public class ExtractorRunnable implements Runnable {

    private final VcfIterator iterator;
    private final SimpleFileWriter writer;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;

    public ExtractorRunnable(
            VcfIterator iterator,
            ChildToParentMap childToParentMap,
            SimpleFileWriter writer
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.writer = writer;

    }

    @Override
    public void run() {

        VcfLine vcfLine;
        while ((vcfLine = iterator.next()) != null) {

            vcfLine.parse();

            String genotypes = processVcfLine(vcfLine);

            String newLine = String.join(
                    "\t",
                    vcfLine.getVariantDescription(),
                    genotypes
            );

            writer.writeLine(newLine);

        }
    }

    private String processVcfLine(VcfLine vcfLine) {

        return childToParentMap.children.stream()
                .parallel()
                .map(
                        childId -> getTransmittedAllele(
                                vcfLine,
                                childId
                        )
                )
                .collect(
                        Collectors.joining("\t")
                );

    }

    private String getTransmittedAllele(VcfLine vcfLine, String childId) {

        String motherId = childToParentMap.getMother(childId);
        String fatherId = childToParentMap.getFather(childId);
        
        int genotypeKid = vcfLine.getGenotype(childId);
        int genotypeMother = vcfLine.getGenotype(motherId);
        int genotypeFather = vcfLine.getGenotype(fatherId);
        
        return Integer.toString(genotypeKid + genotypeMother + genotypeFather);
        
    }

}
