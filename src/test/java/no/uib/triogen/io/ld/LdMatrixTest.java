package no.uib.triogen.io.ld;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.zip.Deflater;
import junit.framework.Assert;
import junit.framework.TestCase;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.model.trio_genotypes.VariantIndex;

/**
 * This class tests the ld matrix reading and writing.
 *
 * @author Marc Vaudel
 */
public class LdMatrixTest extends TestCase {

    public void testParsing() {

        boolean success = false;

        try {

            VariantIndex variantIndex = new VariantIndex();

            File matrixFile = new File("src/test/resources/ld/test.tld");
            LdMatrixWriter writer = new LdMatrixWriter(variantIndex, matrixFile);

            // Create a dummy LD map and store it
            TreeMap<String, ArrayList<R2>> ldMap = new TreeMap<>();

            for (int i = 0; i < 100; i++) {

                ArrayList<R2> r2s = new ArrayList<>(100);

                for (int j = 0; j < 100; j++) {

                    String variantB = "variantB_" + j;
                    double r2Value = ((double) i) / 200 + ((double) j) / 200;

                    variantIndex.add(variantB);

                    int variantBI = variantIndex.getIndex(variantB);

                    R2 r2 = new R2(variantBI, (short) (100 - j), (short) j, (float) r2Value);

                    r2s.add(r2);

                }

                String variant = "variantA_" + i;
                ldMap.put(variant, r2s);

                int variantAI = variantIndex.getIndex(variant);

                writer.addVariant(variantAI, r2s);

            }

            writer.close();

            // Check data retrieval
            LdMatrixReader ldMatrixReader = new LdMatrixReader(matrixFile);

            ArrayList<R2> dummyMapping = ldMatrixReader.getR2("DUMMY");
            Assert.assertTrue(dummyMapping == null);

            for (String variantA : ldMap.keySet()) {

                ArrayList<R2> groundTruth = ldMap.get(variantA);

                ArrayList<R2> r2s = ldMatrixReader.getR2(variantA);
                Assert.assertTrue(r2s != null);
                Assert.assertTrue(r2s.size() == groundTruth.size());

                for (int i = 0; i < r2s.size(); i++) {

                    R2 fileR2 = r2s.get(i);
                    R2 groundTruthR2 = groundTruth.get(i);

                    Assert.assertTrue(fileR2.variantB == groundTruthR2.variantB);
                    Assert.assertTrue(fileR2.alleleA == groundTruthR2.alleleA);
                    Assert.assertTrue(fileR2.alleleB == groundTruthR2.alleleB);

                    Assert.assertTrue(Math.abs(fileR2.r2Value - groundTruthR2.r2Value) <= 1e-6);

                }
            }

            success = true;

        } catch (Throwable t) {

            t.printStackTrace();

        }

        Assert.assertTrue(success);

    }

}
