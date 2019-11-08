package no.uib.triogen.io.vcf;

import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.io.genotypes.vcf.custom.VcfLine;
import java.io.File;
import junit.framework.Assert;
import junit.framework.TestCase;

/**
 *
 * @author Marc Vaudel
 */
public class VcfParsingTest extends TestCase {

    public void testParsing() {
        
        File vcfFile = new File("src/test/resources/vcf/example/test.vcf");

        CustomVcfIterator iterator = new CustomVcfIterator(vcfFile);

        int nLines = 0;

        VcfLine vcfLine;
        while ((vcfLine = iterator.next()) != null) {

            vcfLine.parse();

            String variantDescription = vcfLine.getVariantDescription();
            Assert.assertTrue(variantDescription.equals("1	123	rs123	A	B	12	bla	bla"));
            
            int genotype = vcfLine.getGenotype("SAMPLE1");
            Assert.assertTrue(genotype == 0);
            
            genotype = vcfLine.getGenotype("SAMPLE2");
            Assert.assertTrue(genotype == 2);

            nLines++;

        }

        Assert.assertTrue(nLines == 1);

    }

}
