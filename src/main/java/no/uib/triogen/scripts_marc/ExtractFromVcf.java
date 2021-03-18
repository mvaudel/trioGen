package no.uib.triogen.scripts_marc;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.vcf.iterators.VcfIterator;
import static no.uib.triogen.io.genotypes.vcf.iterators.VcfIterator.getVcfIndexFile;
import no.uib.triogen.io.genotypes.vcf.reader.VcfVariant;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Extract data from vcf files.
 *
 * @author Marc Vaudel
 */
public class ExtractFromVcf {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String vcfFilePath = "/mnt/archive/utils/dbsnp/human_9606_b151_GRCh37p13/All_20180423.vcf.gz";

        System.out.println(Instant.now() + "    Querying " + vcfFilePath + ".");

        String targetChr = "2";
        int bp = 46340137;

        Instant begin = Instant.now();

        File vcfFile = new File(vcfFilePath);
        File indexFile = getVcfIndexFile(vcfFile);
        VCFFileReader reader = new VCFFileReader(vcfFile, indexFile);

        try (CloseableIterator<VariantContext> iterator = reader.query(targetChr, bp, bp)) {

            VariantContext variantContext;
            while ((variantContext = iterator.next()) != null) {

                String alleles = variantContext.getAlleles().stream()
                        .map(
                                allele -> allele.getBaseString()
                        )
                        .collect(
                                Collectors.joining(",")
                        );

                System.out.println("chr: " + variantContext.getContig() + "; id: " + variantContext.getID() + "; id: " + alleles);

            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Done (" + durationSeconds + " s)");

    }
}
