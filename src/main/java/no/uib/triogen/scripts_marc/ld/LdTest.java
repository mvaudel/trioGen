package no.uib.triogen.scripts_marc.ld;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.genotypes.InheritanceUtils;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.processing.ld.LdMatrixComputer;

/**
 * Extracts specific lines of a file.
 *
 * @author Marc Vaudel
 */
public class LdTest {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File genotypesFile = new File("tmp/22.phased.bgen");
            File trioFile = new File("tmp/trio");
            String destinationFile = "tmp/22_test";
            String chromosome = "20";

        ChildToParentMap childToParentMap = ChildToParentMap.fromFile(trioFile);

        HashMap<Integer, char[]> inheritanceMap = InheritanceUtils.getDefaultInheritanceMap(chromosome);
        int defaultMotherPlooidy = InheritanceUtils.getDefaultMotherPloidy(chromosome);
        int defaultFatherPlooidy = InheritanceUtils.getDefaultFatherPloidy(chromosome);

        File logFile = new File("tmp/20_test.log.gz");
        SimpleCliLogger logger = new SimpleCliLogger(logFile, null);
        logger.writeComment("Software", "TrioGen");
        logger.writeComment("Version", TrioGen.getVersion());
        logger.writeComment("Command", "LinearModel");
        logger.writeComment("Arguments", "test");
        logger.writeHeaders();

        LdMatrixComputer computer = new LdMatrixComputer(
                genotypesFile,
                inheritanceMap,
                defaultMotherPlooidy,
                defaultFatherPlooidy,
                childToParentMap,
                destinationFile,
                500000,
                1e-6,
                1e-3,
                2,
                logger
        );

        try {

            computer.run(
                    365
            );

        } catch (Throwable e) {

            e.printStackTrace();

        }

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
