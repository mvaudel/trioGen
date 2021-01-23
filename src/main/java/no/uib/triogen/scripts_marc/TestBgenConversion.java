package no.uib.triogen.scripts_marc;

import java.io.File;
import no.uib.triogen.io.conversion.VcfToBgenConverter;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 *
 *
 * @author Marc Vaudel
 */
public class TestBgenConversion {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File vcfFile = new File("C:\\Github\\trioGen\\tmp\\23.vcf.gz");
        File bgenFile = new File("C:\\Github\\trioGen\\tmp\\23.phased.bgen");

        try {

            VcfToBgenConverter vcfToBgenConverter = new VcfToBgenConverter();
            
            vcfToBgenConverter.convert(vcfFile, bgenFile);

        } catch (Exception e) {

            e.printStackTrace();

        }

    }

}
