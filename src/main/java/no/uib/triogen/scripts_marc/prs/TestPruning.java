package no.uib.triogen.scripts_marc.prs;

import java.io.File;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.processing.prs.PrsPruner;

/**
 * This class runs WLM on meta results. Based on work by RN Beaumont.
 *
 * @author Marc Vaudel
 */
public class TestPruning {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        try {

            String ldMatrixFilePath = "tmp/{contig}_200k.tld";
            File wlmFile = new File("tmp/prs_wlm.gz");
            File destinationFile = new File("tmp/test_prs.gz");
            File logFile = new File("tmp/test_prs_log.gz");

            SimpleCliLogger logger = new SimpleCliLogger(logFile, null);

            PrsPruner prsTrainer = new PrsPruner(
                    wlmFile,
                    ldMatrixFilePath,
                    destinationFile,
                    "moba_id",
                    "chr",
                    "pos",
                    "other_allele",
                    "tested_allele",
                    "beta_wlm_{variable}",
                    "se_wlm_{variable}",
                    "p_wlm_{variable}",
                    new String[]{"child", "mother", "father"},
                    5,
                    0.05,
                    0.8,
                    0.05,
                    logger
            );

            prsTrainer.run();

        } catch (Throwable t) {

            t.printStackTrace();

        }

    }
}
