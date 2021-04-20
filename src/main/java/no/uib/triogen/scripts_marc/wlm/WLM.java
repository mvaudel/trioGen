package no.uib.triogen.scripts_marc.wlm;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This class runs WLM on meta results. Based on work by RN Beaumont.
 *
 * @author Marc Vaudel
 */
public class WLM {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        double childMotherOverlap = 0.2653;
        double childFatherOverlap = 0.197;
        double motherFatherOverlap = 0.0082;

        File gwasChild = new File("/mnt/work/marc/moba/pwbw/prs/meta/child/child_prs1_rsid_clean.tbl.gz");
        File gwasMother = new File("/mnt/work/marc/moba/pwbw/prs/meta/mother/mother_prs1_rsid_clean.tbl.gz");
        File gwasFather = new File("/mnt/work/marc/moba/pwbw/prs/meta/father/father_prs1_rsid_clean.tbl.gz");

        File resultFile = new File("/mnt/work/marc/moba/pwbw/prs/wlm/prs_wlm.gz");

        HashMap<String, String[]> variantInfoMap = new HashMap<>();
        HashMap<String, String[]> allelesMap = new HashMap<>();
        HashMap<String, double[]> frequencyMap = new HashMap<>();
        HashMap<String, double[]> betaMap = new HashMap<>();
        HashMap<String, double[]> seMap = new HashMap<>();
        HashMap<String, double[]> pMap = new HashMap<>();
        HashMap<String, double[]> nMap = new HashMap<>();

        // Load child

        System.out.println(Instant.now() + "    Loading child results.");
        
        try (SimpleFileReader reader = SimpleFileReader.getFileReader(gwasChild)) {

            String line = reader.readLine();

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                String id = lineSplit[0];
                String rsid = lineSplit[1];
                String chr = lineSplit[2];
                String pos = lineSplit[3];
                String testedAllele = lineSplit[4];
                String otherAllele = lineSplit[5];
                double testedAlleleFrequency = Double.parseDouble(lineSplit[6]);
                double beta = Double.parseDouble(lineSplit[10]);
                double se = Double.parseDouble(lineSplit[11]);
                double p = Double.parseDouble(lineSplit[12]);
                double n = Double.parseDouble(lineSplit[18]);

                if (beta < 0) {

                    beta = -beta;

                    String temp = testedAllele;
                    testedAllele = otherAllele;
                    otherAllele = temp;

                    testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                }

                variantInfoMap.put(id, new String[]{rsid, chr, pos});
                allelesMap.put(id, new String[]{otherAllele, testedAllele});
                frequencyMap.put(id, new double[]{testedAlleleFrequency, Double.NaN, Double.NaN});
                betaMap.put(id, new double[]{beta, Double.NaN, Double.NaN});
                seMap.put(id, new double[]{se, Double.NaN, Double.NaN});
                pMap.put(id, new double[]{p, Double.NaN, Double.NaN});
                nMap.put(id, new double[]{n, 0, 0});

            }
        }

        // Load mother

        System.out.println(Instant.now() + "    Loading mother results.");
        
        try (SimpleFileReader reader = SimpleFileReader.getFileReader(gwasMother)) {

            String line = reader.readLine();

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                String id = lineSplit[0];

                if (variantInfoMap.containsKey(id)) {

                    String testedAllele = lineSplit[4];
                    String otherAllele = lineSplit[5];
                    double testedAlleleFrequency = Double.parseDouble(lineSplit[6]);
                    double beta = Double.parseDouble(lineSplit[10]);
                    double se = Double.parseDouble(lineSplit[11]);
                    double p = Double.parseDouble(lineSplit[12]);
                    double n = Double.parseDouble(lineSplit[18]);

                    if (beta < 0) {

                        beta = -beta;

                        String temp = testedAllele;
                        testedAllele = otherAllele;
                        otherAllele = temp;

                        testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                    }

                    String[] alleles = allelesMap.get(id);
                    double[] betas = betaMap.get(id);
                    double[] ses = seMap.get(id);
                    double[] testedAlleleFrequencies = frequencyMap.get(id);
                    double[] ps = pMap.get(id);
                    double[] ns = nMap.get(id);

                    if (ps[0] <= p) {

                        if (alleles[0].equals(testedAllele) && alleles[1].equals(otherAllele)) {

                            beta = -beta;
                            testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                        } else if (!alleles[0].equals(otherAllele) || !alleles[1].equals(testedAllele)) {

                            throw new IllegalArgumentException("Allele mismatch for variant " + id + ".");

                        }
                    } else {

                        if (alleles[0].equals(testedAllele) && alleles[1].equals(otherAllele)) {

                            betas[0] = -betas[0];
                            testedAlleleFrequencies[0] = 1.0 - testedAlleleFrequencies[0];

                        } else if (!alleles[0].equals(otherAllele) || !alleles[1].equals(testedAllele)) {

                            throw new IllegalArgumentException("Allele mismatch for variant " + id + ".");

                        }
                    }

                    betas[1] = beta;
                    ses[1] = se;
                    ps[1] = p;
                    testedAlleleFrequencies[1] = testedAlleleFrequency;
                    ns[1] = n;

                }
            }
        }

        // Load father

        System.out.println(Instant.now() + "    Loading father results.");
        
        try (SimpleFileReader reader = SimpleFileReader.getFileReader(gwasFather)) {

            String line = reader.readLine();

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                String id = lineSplit[0];

                if (variantInfoMap.containsKey(id)) {

                    String testedAllele = lineSplit[4];
                    String otherAllele = lineSplit[5];
                    double testedAlleleFrequency = Double.parseDouble(lineSplit[6]);
                    double beta = Double.parseDouble(lineSplit[10]);
                    double se = Double.parseDouble(lineSplit[11]);
                    double p = Double.parseDouble(lineSplit[12]);
                    double n = Double.parseDouble(lineSplit[18]);

                    if (beta < 0) {

                        beta = -beta;

                        String temp = testedAllele;
                        testedAllele = otherAllele;
                        otherAllele = temp;

                        testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                    }

                    String[] alleles = allelesMap.get(id);
                    double[] betas = betaMap.get(id);
                    double[] ses = seMap.get(id);
                    double[] testedAlleleFrequencies = frequencyMap.get(id);
                    double[] ps = pMap.get(id);
                    double[] ns = nMap.get(id);

                    if (ps[0] <= ps[1] && ps[0] <= p) {

                        if (alleles[0].equals(testedAllele) && alleles[1].equals(otherAllele)) {

                            beta = -beta;
                            testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                        } else if (!alleles[0].equals(otherAllele) || !alleles[1].equals(testedAllele)) {

                            throw new IllegalArgumentException("Allele mismatch for variant " + id + ".");

                        }
                    } else if (ps[1] <= ps[0] && ps[1] <= p) {

                        if (alleles[0].equals(testedAllele) && alleles[1].equals(otherAllele)) {

                            beta = -beta;
                            testedAlleleFrequency = 1.0 - testedAlleleFrequency;

                        } else if (!alleles[0].equals(otherAllele) || !alleles[1].equals(testedAllele)) {

                            throw new IllegalArgumentException("Allele mismatch for variant " + id + ".");

                        }
                    } else {

                        if (alleles[0].equals(testedAllele) && alleles[1].equals(otherAllele)) {

                            betas[0] = -betas[0];
                            betas[1] = -betas[1];
                            testedAlleleFrequencies[0] = 1.0 - testedAlleleFrequencies[0];
                            testedAlleleFrequencies[1] = 1.0 - testedAlleleFrequencies[1];

                        } else if (!alleles[0].equals(otherAllele) || !alleles[1].equals(testedAllele)) {

                            throw new IllegalArgumentException("Allele mismatch for variant " + id + ".");

                        }
                    }

                    betas[2] = beta;
                    ses[2] = se;
                    ps[2] = p;
                    testedAlleleFrequencies[2] = testedAlleleFrequency;
                    ns[2] = n;

                }
            }
        }

        // Compute WLM and export to file

        System.out.println(Instant.now() + "    Computing WLM");
        
        try (SimpleFileWriter writer = new SimpleFileWriter(resultFile, true)) {

            writer.writeLine(
                    "snp", "rsid", "chr", "pos", "tested_allele", "other_allele",
                    "tested_allele_freq_child", "tested_allele_freq_mother", "tested_allele_freq_father",
                    "beta_child", "se_child", "p_child",
                    "beta_mother", "se_mother", "p_mother",
                    "beta_father", "se_father", "p_father",
                    "beta_wlm_child", "se_wlm_child",
                    "beta_wlm_mother", "se_wlm_mother",
                    "beta_wlm_father", "se_wlm_father"
            );

            for (String id : variantInfoMap.keySet()) {

                double[] ps = pMap.get(id);

                if (!Double.isNaN(ps[0]) && !Double.isNaN(ps[1]) && !Double.isNaN(ps[2])) {

                    String[] variantInfo = variantInfoMap.get(id);
                    String[] alleles = allelesMap.get(id);
                    double[] testedAlleleFrequencies = frequencyMap.get(id);
                    double[] betas = betaMap.get(id);
                    double[] ses = seMap.get(id);
                    double[] ns = nMap.get(id);

                    double[] ses2 = Arrays.stream(ses)
                            .map(
                                    se -> se * se
                            )
                            .toArray();

                    double wlmBetaChild = 2 * betas[0] - betas[1] - betas[1];
                    double wlmSeChild = Math.sqrt(
                            4 * ses2[0] + ses2[1] + ses2[2]
                            - 4 * childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                            - 4 * childFatherOverlap * Math.sqrt(ses[0] * ses2[2])
                            + 2 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2])
                    );
                    double wlmBetaMother = (3 * betas[1] - 2 * betas[0] + betas[2]) / 2;
                    double wlmSeMother = Math.sqrt(
                            9 * ses2[1] / 4 + ses2[0] + ses2[2] / 4
                            - 3 * childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                            - childFatherOverlap * Math.sqrt(ses[0] * ses2[2])
                            + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
                    );
                    double wlmBetaFather = (3 * betas[2] - 2 * betas[0] + betas[1]) / 2;
                    double wlmSeFather = Math.sqrt(
                            9 * ses2[2] / 4 + ses2[0] + ses2[1] / 4
                            - childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                            - 3 * childFatherOverlap * Math.sqrt(ses[0] * ses2[2])
                            + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
                    );

                    writer.writeLine(
                            id, variantInfo[0], variantInfo[1], variantInfo[2], alleles[1], alleles[0],
                            Double.toString(testedAlleleFrequencies[0]), Double.toString(testedAlleleFrequencies[1]), Double.toString(testedAlleleFrequencies[2]),
                            Double.toString(betas[0]), Double.toString(ses[0]), Double.toString(ps[0]),
                            Double.toString(betas[1]), Double.toString(ses[1]), Double.toString(ps[1]),
                            Double.toString(betas[2]), Double.toString(ses[2]), Double.toString(ps[2]),
                            Double.toString(wlmBetaChild), Double.toString(wlmSeChild),
                            Double.toString(wlmBetaMother), Double.toString(wlmSeMother),
                            Double.toString(wlmBetaFather), Double.toString(wlmSeFather)
                    );

                }
            }
        }
    }
}
