package no.uib.triogen.io.ld;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.zip.Deflater;
import junit.framework.Assert;
import junit.framework.TestCase;
import no.uib.triogen.model.trio_genotypes.VariantIndex;

/**
 * This class tests the ld matrix reading and writing.
 *
 * @author Marc Vaudel
 */
public class LdMatrixTest extends TestCase {

    public void testParsing() {
        
        Deflater deflater = new Deflater(Deflater.BEST_COMPRESSION, true);

        boolean success = false;

        try {

            VariantIndex variantIndex = new VariantIndex();

            File matrixFile = new File("src/test/resources/ld/test.tld");
            LdMatrixWriter writer = new LdMatrixWriter(variantIndex, matrixFile);

            // Create a dummy LD map and store it
            TreeMap<String, HashMap<String, Double>> ldMap = new TreeMap<>();

            for (int i = 0; i < 100; i++) {

                HashMap<String, Double> variantMap = new HashMap<>(100);
                ArrayList<Integer> variantBs = new ArrayList<>(100);
                ArrayList<Double> r2s = new ArrayList<>(100);

                for (int j = 0; j < 100; j++) {

                    String variantB = "variantB_" + j;
                    double r2 = ((double) i) / 200 + ((double) j) / 200;
                    variantMap.put(variantB, r2);

                    variantIndex.add(variantB);

                    int variantBI = variantIndex.getIndex(variantB);
                    variantBs.add(variantBI);
                    r2s.add(r2);

                }

                String variant = "variantA_" + i;
                ldMap.put(variant, variantMap);

                int variantAI = variantIndex.getIndex(variant);

                writer.addVariant(variantAI, variantBs, r2s, deflater);

            }

            writer.close();

            // Check data retrieval
            LdMatrixReader ldMatrixReader = new LdMatrixReader(matrixFile);

            HashMap<String, Double> dummyMapping = ldMatrixReader.getR2("DUMMY");
            Assert.assertTrue(dummyMapping == null);

            for (String variantA : ldMap.keySet()) {

                HashMap<String, Double> groundTruth = ldMap.get(variantA);

                HashMap<String, Double> ldMapping = ldMatrixReader.getR2(variantA);
                Assert.assertTrue(ldMapping != null);

                for (Entry<String, Double> entry : groundTruth.entrySet()) {

                    String variantB = entry.getKey();
                    double r2GT = entry.getValue();

                    Assert.assertTrue(ldMapping.containsKey(variantB));

                    double r2File = ldMapping.get(variantB);

                    Assert.assertTrue(r2GT == r2File);

                    dummyMapping = ldMatrixReader.getR2(variantB);
                    Assert.assertTrue(dummyMapping == null);

                }
            }

            success = true;

        } catch (Throwable t) {

            t.printStackTrace();

        }

        Assert.assertTrue(success);

    }

}
