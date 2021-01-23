package no.uib.triogen.scripts_marc;

import java.io.File;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 *
 *
 * @author Marc Vaudel
 */
public class TestBgen {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File bgenFile = new File("C:\\Github\\trioGen\\tmp\\23.phased.bgen");

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(new File("C:\\Github\\trioGen\\tmp\\trio"));

        try {

            BgenIndex index = BgenIndex.getBgenIndex(bgenFile);

            BgenFileReader reader = new BgenFileReader(bgenFile, index, null, 0);

            int phased = 0;

            for (int i = 0; i < reader.getNVariants(); i++) {

                BgenVariantData variantData = reader.getVariantData(i);

                try {

                    variantData.parse(childToParentMap);
                    phased++;

                } catch (Exception e) {

                }

            }

            System.out.println(phased + " variants phased.");

        } catch (Exception e) {

            e.printStackTrace();

        }

    }

}
