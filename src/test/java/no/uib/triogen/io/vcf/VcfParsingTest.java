package no.uib.triogen.io.vcf;

import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.io.genotypes.vcf.custom.VcfLine;
import java.io.File;
import java.util.HashMap;
import junit.framework.Assert;
import junit.framework.TestCase;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 *
 * @author Marc Vaudel
 */
public class VcfParsingTest extends TestCase {

    public void testParsing() {

        String[] children = new String[]{"SAMPLE1"};
        HashMap<String, String> childToMotherMap = new HashMap<>();
        childToMotherMap.put("SAMPLE1", "SAMPLE2");
        HashMap<String, String> childToFatherMap = new HashMap<>();
        childToFatherMap.put("SAMPLE1", "SAMPLE3");
        
        ChildToParentMap childToParentMap = new ChildToParentMap(children, childToMotherMap, childToFatherMap);
        
        File vcfFile = new File("src/test/resources/vcf/example/test.vcf");

        CustomVcfIterator iterator = new CustomVcfIterator(
                vcfFile
        );

        int nLines = 0;

        VcfLine vcfLine;
        while ((vcfLine = iterator.next()) != null) {

            if (nLines == 0) {

                vcfLine.parse(childToParentMap);

                String variantId = vcfLine.getVariantID();
                Assert.assertTrue(variantId.equals("rs123"));

                String contig = vcfLine.getContig();
                Assert.assertTrue(contig.equals("1"));

                int bp = vcfLine.getBp();
                Assert.assertTrue(bp == 123);

                String ref = vcfLine.getRef();
                Assert.assertTrue(ref.equals("A"));

                String alt = vcfLine.getAlt();
                Assert.assertTrue(alt.equals("C"));

                boolean genotyped = vcfLine.genotyped();
                Assert.assertTrue(!genotyped);

                int genotype = vcfLine.getGenotype("SAMPLE1");
                Assert.assertTrue(genotype == 0);

                float[] dosages = vcfLine.getDosages("SAMPLE1");
                double[] groundTruth = new double[]{0.9025, 0.095, 0.0025};
                Assert.assertTrue(dosages.length == 3);

                for (int i = 0; i < 3; i++) {

                    double expected = groundTruth[0];
                    double found = dosages[0];
                    Assert.assertTrue(Math.abs(expected - found) < 1e-5);

                }

                genotype = vcfLine.getGenotype("SAMPLE2");
                Assert.assertTrue(genotype == 2);

                dosages = vcfLine.getDosages("SAMPLE2");
                groundTruth = new double[]{0.1, 0.8, 0.1};
                Assert.assertTrue(dosages.length == 3);

                for (int i = 0; i < 3; i++) {

                    double expected = groundTruth[0];
                    double found = dosages[0];
                    Assert.assertTrue(Math.abs(expected - found) < 1e-5);

                }

            } else if (nLines == 1) {

                vcfLine.parse(childToParentMap);

                String variantId = vcfLine.getVariantID();
                Assert.assertTrue(variantId.equals("rs234"));

                String contig = vcfLine.getContig();
                Assert.assertTrue(contig.equals("2"));

                int bp = vcfLine.getBp();
                Assert.assertTrue(bp == 234);

                String ref = vcfLine.getRef();
                Assert.assertTrue(ref.equals("T"));

                String alt = vcfLine.getAlt();
                Assert.assertTrue(alt.equals("G"));

                boolean genotyped = vcfLine.genotyped();
                Assert.assertTrue(genotyped);

                int genotype = vcfLine.getGenotype("SAMPLE1");
                Assert.assertTrue(genotype == 1);

                float[] dosages = vcfLine.getDosages("SAMPLE1");
                double[] groundTruth = new double[]{0.0, 0.0, 1.0};
                Assert.assertTrue(dosages.length == 3);

                for (int i = 0; i < 3; i++) {

                    double expected = groundTruth[0];
                    double found = dosages[0];
                    Assert.assertTrue(Math.abs(expected - found) < 1e-5);

                }

                genotype = vcfLine.getGenotype("SAMPLE2");
                Assert.assertTrue(genotype == 3);

                dosages = vcfLine.getDosages("SAMPLE2");
                groundTruth = new double[]{1.0, 0.0, 0.0};
                Assert.assertTrue(dosages.length == 3);

                for (int i = 0; i < 3; i++) {

                    double expected = groundTruth[0];
                    double found = dosages[0];
                    Assert.assertTrue(Math.abs(expected - found) < 1e-5);

                }
            }

            nLines++;

        }

        Assert.assertTrue(nLines == 2);

    }

}
