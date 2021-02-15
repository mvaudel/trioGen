package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.genotypes.InheritanceUtils;
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

            System.out.println(Instant.now() + " Creating index...");

            BgenIndex index = BgenIndex.getBgenIndex(bgenFile);

            System.out.println(Instant.now() + " Index created.");

            int defaultMotherPlooidy = InheritanceUtils.getDefaultMotherPloidy("23");
            int defaultFatherPlooidy = InheritanceUtils.getDefaultFatherPloidy("23");

            BgenFileReader reader = new BgenFileReader(bgenFile, index, InheritanceUtils.getDefaultInheritanceMap("23"), defaultMotherPlooidy, defaultFatherPlooidy);

            int phased = 0;
            int previousProgress = 0;

            for (int i = 0; i < index.variantInformationArray.length; i++) {

                double progress = (100.0 * i) / index.variantInformationArray.length;

                if (progress >= previousProgress + 1) {

                    System.out.println(Instant.now() + " Parsing variants... " + i + " of " + index.variantInformationArray.length + " (" + ((int) progress) + "%)");

                    previousProgress = (int) progress;

                }

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
