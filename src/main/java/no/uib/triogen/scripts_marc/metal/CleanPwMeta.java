package no.uib.triogen.scripts_marc.metal;

import java.io.File;
import java.time.Instant;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * Filters and arranges the results of the meta for Pw.
 *
 * @author Marc Vaudel
 */
public class CleanPwMeta {

    public static final double MIN_SAMPLE_SIZE = 2500;
    public static final double MIN_HET_P = 5e-8;
    public static final double MIN_STUDIES = 2;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            String rawFilePath = "/mnt/work/marc/moba/pwbw/prs_no_moba/meta/{geno}/{geno}_prs_no_moba_1_rsid.tbl.gz";
            String cleanFilePath = "/mnt/work/marc/moba/pwbw/prs_no_moba/meta/{geno}/{geno}_prs_no_moba_1_rsid_clean.tbl.gz";

            String[] genos = new String[]{"child", "mother", "father"};

            for (String geno : genos) {

                File rawFile = new File(rawFilePath.replace("{geno}", geno));
                File cleanFile = new File(cleanFilePath.replace("{geno}", geno));
                
                System.out.println(Instant.now() + "    Processing " + rawFile + ".");

                try (SimpleFileReader reader = SimpleFileReader.getFileReader(rawFile)) {

                    try (SimpleFileWriter writer = new SimpleFileWriter(cleanFile, true)) {

                        writer.writeLine(
                                "snp",
                                "rsid",
                                "chr",
                                "pos",
                                "tested_allele",
                                "other_allele",
                                "tested_allele_freq",
                                "tested_allele_freq_se",
                                "tested_allele_freq_min",
                                "tested_allele_freq_max",
                                "beta",
                                "se",
                                "p",
                                "direction",
                                "het_i2",
                                "het_chi2",
                                "het_df",
                                "het_p",
                                "total_sample_size"
                        );

                        String line = reader.readLine();

                        while ((line = reader.readLine()) != null) {

                            String[] lineSplit = line.split("\t");

                            String snp = lineSplit[0];
                            String rsId = lineSplit[1];
                            String testedAllele = lineSplit[2].toUpperCase();
                            String otherAllele = lineSplit[3].toUpperCase();
                            String af = lineSplit[4];
                            String afSe = lineSplit[5];
                            String afMin = lineSplit[6];
                            String afMax = lineSplit[7];
                            String beta = lineSplit[8];
                            String se = lineSplit[9];
                            String p = lineSplit[10];
                            String direction = lineSplit[11];
                            String hetI2 = lineSplit[12];
                            String hetChi2 = lineSplit[13];
                            String hetDf = lineSplit[14];
                            String hetP = lineSplit[15];
                            String n = lineSplit[16];
                            
                            String[] tmpSplit = snp.split(":");
                            String chr = tmpSplit[0];
                            
                            tmpSplit = tmpSplit[1].split("_");
                            String pos = tmpSplit[0];

                            int nStudies = 0;

                            for (int i = 0; i < direction.length(); i++) {

                                if (direction.charAt(i) != '?') {

                                    nStudies++;

                                }
                            }

                            if (Double.parseDouble(n) > MIN_SAMPLE_SIZE
                                    && nStudies >= MIN_STUDIES
                                    && Double.parseDouble(hetP) >= MIN_HET_P) {

                                writer.writeLine(
                                        snp,
                                        rsId,
                                        chr,
                                        pos,
                                        testedAllele,
                                        otherAllele,
                                        af,
                                        afSe,
                                        afMin,
                                        afMax,
                                        beta,
                                        se,
                                        p,
                                        direction,
                                        hetI2,
                                        hetChi2,
                                        hetDf,
                                        hetP,
                                        n
                                );

                            } else {
                                
                                System.out.println(snp + " " + n + "," + nStudies + "," + hetP);
                                
                            }
                        }
                    }
                }
            }

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
