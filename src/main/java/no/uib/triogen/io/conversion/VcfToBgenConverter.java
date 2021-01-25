package no.uib.triogen.io.conversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.writer.BgenWriter;
import no.uib.triogen.io.genotypes.vcf.reader.VcfIterator;
import no.uib.triogen.io.genotypes.vcf.reader.VcfVariant;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * This class converts a vcf file to a bgen file.
 *
 * @author Marc Vaudel
 */
public class VcfToBgenConverter {

    public void convert(
            File vcfFile,
            File bgenFile
    ) throws IOException {

        try (VcfIterator vcfIterator = new VcfIterator(vcfFile)) {

            try (BgenWriter bgenWriter = new BgenWriter(bgenFile, BgenIndex.getDefaultIndexFile(bgenFile))) {

                bgenWriter.initiate(vcfIterator.getSamplesIds());

                VcfVariant vcfVariant;

                while ((vcfVariant = vcfIterator.next()) != null) {

                    VariantInformation variantInformation = vcfVariant.getVariantInformation();
                    ArrayList<String[]> genotypes = vcfVariant.getGenotypes(vcfIterator.getSamplesIds());

                    bgenWriter.addVariant(variantInformation, genotypes);

                }

                bgenWriter.finalize();

            }
        }
    }
}
