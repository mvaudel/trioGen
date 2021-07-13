package no.uib.triogen.scripts_marc.wlm;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.utils.Utils;

/**
 * This class runs WLM on meta results. Based on work by RN Beaumont.
 *
 * @author Marc Vaudel
 */
public class WLM {

    private static double childMotherOverlap = 0.2653;
    private static double childFatherOverlap = 0.197;
    private static double motherFatherOverlap = 0.0082;

    private static int currentProgress = -1;
    private static int totalProgress = 0;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File gwasChild = new File("/mnt/work/marc/moba/pwbw/prs_no_moba/meta/child/child_prs_no_moba_1_rsid_clean.tbl.gz");
        File gwasMother = new File("/mnt/work/marc/moba/pwbw/prs_no_moba/meta/mother/mother_prs_no_moba_1_rsid_clean.tbl.gz");
        File gwasFather = new File("/mnt/work/marc/moba/pwbw/prs_no_moba/meta/father/father_prs_no_moba_1_rsid_clean.tbl.gz");
        File resultFile = new File("/mnt/work/marc/moba/pwbw/prs_no_moba/wlm/prs_wlm.gz");

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
                            allelesMap.put(id, new String[]{otherAllele, testedAllele});

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
                            allelesMap.put(id, new String[]{otherAllele, testedAllele});

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
                    "snp", "moba_id", "rsid", "chr", "pos", "tested_allele", "other_allele",
                    "tested_allele_freq_child", "tested_allele_freq_mother", "tested_allele_freq_father",
                    "beta_child", "se_child", "p_child", "n_child",
                    "beta_mother", "se_mother", "p_mother", "n_mother",
                    "beta_father", "se_father", "p_father", "n_father",
                    "beta_wlm_child", "se_wlm_child", "p_wlm_child",
                    "beta_wlm_mother", "se_wlm_mother", "p_wlm_mother",
                    "beta_wlm_father", "se_wlm_father", "p_wlm_father"
            );

            ArrayList<String> variantIds = new ArrayList<>(variantInfoMap.keySet());
            totalProgress = variantIds.size();
            
            long start = Instant.now().getEpochSecond();

            variantIds.stream()
                    .parallel()
                    .forEach(
                            id -> processVariant(
                                    id,
                                    variantInfoMap,
                                    allelesMap,
                                    frequencyMap,
                                    betaMap,
                                    seMap,
                                    pMap,
                                    nMap,
                                    writer,
                                    start
                            )
                    );
        }
    }

    private static void processVariant(
            String id,
            HashMap<String, String[]> variantInfoMap,
            HashMap<String, String[]> allelesMap,
            HashMap<String, double[]> frequencyMap,
            HashMap<String, double[]> betaMap,
            HashMap<String, double[]> seMap,
            HashMap<String, double[]> pMap,
            HashMap<String, double[]> nMap,
            SimpleFileWriter writer,
            long start
    ) {

        currentProgress++;

        if (currentProgress % 1000 == 0) {

            double progress = currentProgress / totalProgress;
            
            long elapsedTimeSeconds = Instant.now().getEpochSecond() - start;
            
            double remainingTimeSeconds = Math.round(((double) totalProgress) * elapsedTimeSeconds / currentProgress - elapsedTimeSeconds);
            
            double progressDisplay = (10.0 * Math.round(progress) / 10);
            
            System.out.println(Instant.now() + "    Computing WLM - " + currentProgress + " of " + totalProgress + " (" + progressDisplay + " %, ETA: " + remainingTimeSeconds + " s)");

        }

        double[] betas = betaMap.get(id);
        double[] ses = seMap.get(id);

        if (!Double.isNaN(betas[0]) && !Double.isNaN(betas[1]) && !Double.isNaN(betas[2])
                && !Double.isNaN(ses[0]) && !Double.isNaN(ses[1]) && !Double.isNaN(ses[2])) {

            String[] variantInfo = variantInfoMap.get(id);
            String[] alleles = allelesMap.get(id);
            double[] testedAlleleFrequencies = frequencyMap.get(id);
            double[] ps = pMap.get(id);
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
                    - 4 * childFatherOverlap * Math.sqrt(ses2[0] * ses2[2])
                    + 2 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2])
            );
            double wlmPChild = Utils.computeBetaSignificance(wlmBetaChild, wlmSeChild);
            
            double wlmBetaMother = (3 * betas[1] - 2 * betas[0] + betas[2]) / 2;
            double wlmSeMother = Math.sqrt(
                    9 * ses2[1] / 4 + ses2[0] + ses2[2] / 4
                    - 3 * childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                    - childFatherOverlap * Math.sqrt(ses2[0] * ses2[2])
                    + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
            );
            double wlmPMother = Utils.computeBetaSignificance(wlmBetaMother, wlmSeMother);
            
            double wlmBetaFather = (3 * betas[2] - 2 * betas[0] + betas[1]) / 2;
            double wlmSeFather = Math.sqrt(
                    9 * ses2[2] / 4 + ses2[0] + ses2[1] / 4
                    - childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                    - 3 * childFatherOverlap * Math.sqrt(ses2[0] * ses2[2])
                    + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
            );
            double wlmPFather = Utils.computeBetaSignificance(wlmBetaFather, wlmSeFather);
            
            String[] sortedAlleles = Arrays.copyOf(alleles, alleles.length);
            Arrays.sort(sortedAlleles);
            
            String mobaId = String.join("_", variantInfo[1], variantInfo[2], sortedAlleles[0], sortedAlleles[1]);
            
            writer.writeLine(
                    id, mobaId, variantInfo[0], variantInfo[1], variantInfo[2], alleles[1], alleles[0],
                    Double.toString(testedAlleleFrequencies[0]), Double.toString(testedAlleleFrequencies[1]), Double.toString(testedAlleleFrequencies[2]),
                    Double.toString(betas[0]), Double.toString(ses[0]), Double.toString(ps[0]), Double.toString(ns[0]),
                    Double.toString(betas[1]), Double.toString(ses[1]), Double.toString(ps[1]), Double.toString(ns[1]),
                    Double.toString(betas[2]), Double.toString(ses[2]), Double.toString(ps[2]), Double.toString(ns[2]),
                    Double.toString(wlmBetaChild), Double.toString(wlmSeChild), Double.toString(wlmPChild),
                    Double.toString(wlmBetaMother), Double.toString(wlmSeMother), Double.toString(wlmPMother),
                    Double.toString(wlmBetaFather), Double.toString(wlmSeFather), Double.toString(wlmPFather)
            );

        }
    }
}
