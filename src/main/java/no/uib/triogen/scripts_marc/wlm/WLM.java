package no.uib.triogen.scripts_marc.wlm;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.stream.IntStream;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import org.apache.commons.math3.distribution.FDistribution;

/**
 * This class runs WLM on meta results. Based on work by RN Beaumont.
 *
 * @author Marc Vaudel
 */
public class WLM {

    private static double childMotherOverlap = 0.2653;
    private static double childFatherOverlap = 0.197;
    private static double motherFatherOverlap = 0.0082;
    
    private static final int nTriosSimulation = 100000;
    private static final int nTriosPValue = 50000;
    
    private static int lastprogress = -1;
    private static int currentProgress = -1;
    private static int totalProgress = 0;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File gwasChild = new File("/mnt/archive/marc/moba/pwbw/prs/meta/child/child_prs1_rsid_clean.tbl.gz");
        File gwasMother = new File("/mnt/archive/marc/moba/pwbw/prs/meta/mother/mother_prs1_rsid_clean.tbl.gz");
        File gwasFather = new File("/mnt/archive/marc/moba/pwbw/prs/meta/father/father_prs1_rsid_clean.tbl.gz");
        File resultFile = new File("/mnt/archive/marc/moba/pwbw/prs/meta/father/prs_wlm_2.gz");

//        File gwasChild = new File("C:\\Projects\\placenta_weight\\meta_results\\child_prs1_rsid_clean.tbl.gz");
//        File gwasMother = new File("C:\\Projects\\placenta_weight\\meta_results\\mother_prs1_rsid_clean.tbl.gz");
//        File gwasFather = new File("C:\\Projects\\placenta_weight\\meta_results\\father_prs1_rsid_clean.tbl.gz");
//        File resultFile = new File("C:\\Projects\\placenta_weight\\meta_results\\prs_wlm_2.gz");

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
                    "snp", "rsid", "chr", "pos", "tested_allele", "other_allele",
                    "tested_allele_freq_child", "tested_allele_freq_mother", "tested_allele_freq_father",
                    "beta_child", "se_child", "p_child", "n_child",
                    "beta_mother", "se_mother", "p_mother", "n_mother",
                    "beta_father", "se_father", "p_father", "n_father",
                    "beta_wlm_child", "se_wlm_child",
                    "beta_wlm_mother", "se_wlm_mother",
                    "beta_wlm_father", "se_wlm_father",
                    "variance_explained_50", "variance_explained_5", "variance_explained_95", "variance_explained_n",
                    "p_50", "p_5", "p_95", "p_n"
            );

            ArrayList<String> variantIds = new ArrayList<>(variantInfoMap.keySet());
            totalProgress = variantIds.size();

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
                                    writer
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
            SimpleFileWriter writer
    ) {
        
        double progress = 1000.0 * currentProgress / totalProgress;
        
        currentProgress++;
        
        if (progress > lastprogress + 1) {
            
            lastprogress = (int) progress;
            
            double progressDisplay = ((double) lastprogress)/10;
            
            System.out.println(Instant.now() + "    Computing WLM - progress " + progressDisplay + " %");
            
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
            double wlmBetaMother = (3 * betas[1] - 2 * betas[0] + betas[2]) / 2;
            double wlmSeMother = Math.sqrt(
                    9 * ses2[1] / 4 + ses2[0] + ses2[2] / 4
                    - 3 * childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                    - childFatherOverlap * Math.sqrt(ses2[0] * ses2[2])
                    + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
            );
            double wlmBetaFather = (3 * betas[2] - 2 * betas[0] + betas[1]) / 2;
            double wlmSeFather = Math.sqrt(
                    9 * ses2[2] / 4 + ses2[0] + ses2[1] / 4
                    - childMotherOverlap * Math.sqrt(ses2[0] * ses2[1])
                    - 3 * childFatherOverlap * Math.sqrt(ses2[0] * ses2[2])
                    + 3 * motherFatherOverlap * Math.sqrt(ses2[1] * ses2[2]) / 2
            );
            
            double[] trioStats = getTrioSummaryStatisticsRepeated(
                    wlmBetaChild, 
                    wlmBetaMother, 
                    wlmBetaFather, 
                    wlmSeChild, 
                    wlmSeMother, 
                    wlmSeFather, 
                    testedAlleleFrequencies[0], 
                    testedAlleleFrequencies[1], 
                    testedAlleleFrequencies[2]
            );

            writer.writeLine(
                    id, variantInfo[0], variantInfo[1], variantInfo[2], alleles[1], alleles[0],
                    Double.toString(testedAlleleFrequencies[0]), Double.toString(testedAlleleFrequencies[1]), Double.toString(testedAlleleFrequencies[2]),
                    Double.toString(betas[0]), Double.toString(ses[0]), Double.toString(ps[0]), Double.toString(ns[0]),
                    Double.toString(betas[1]), Double.toString(ses[1]), Double.toString(ps[1]), Double.toString(ns[1]),
                    Double.toString(betas[2]), Double.toString(ses[2]), Double.toString(ps[2]), Double.toString(ns[2]),
                    Double.toString(wlmBetaChild), Double.toString(wlmSeChild),
                    Double.toString(wlmBetaMother), Double.toString(wlmSeMother),
                    Double.toString(wlmBetaFather), Double.toString(wlmSeFather),
                    Double.toString(trioStats[0]), Double.toString(trioStats[1]), Double.toString(trioStats[2]), Integer.toString(nTriosSimulation),
                    Double.toString(trioStats[3]), Double.toString(trioStats[4]), Double.toString(trioStats[5]), Integer.toString(nTriosPValue)
            );

        }
    }

    private static double[] getTrioSummaryStatisticsRepeated(
            double childBeta,
            double motherBeta,
            double fatherBeta,
            double childSe,
            double motherSe,
            double fatherSe,
            double mafChild,
            double mafMother,
            double mafFather
    ) {
        
        ArrayList<Double> explainedVariance = new ArrayList<>(101);
        ArrayList<Double> p = new ArrayList<>(101);
        
        for (int i = 0 ; i < 101 ; i++) {
            
            double[] runResults = getTrioSummaryStatistics(
                    childBeta, 
                    motherBeta, 
                    fatherBeta, 
                    childSe, 
                    motherSe, 
                    fatherSe, 
                    mafChild, 
                    mafMother, 
                    mafFather
            );
            
            explainedVariance.add(runResults[0]);
            p.add(runResults[1]);
            
        }
        
        Collections.sort(explainedVariance);
        Collections.sort(p);
        
        double var5 = explainedVariance.get(0);
        double var50 = explainedVariance.get(50);
        double var95 = explainedVariance.get(100);
        
        
        double p5 = p.get(0);
        double p50 = p.get(50);
        double p95 = p.get(100);
        
        return new double[]{var50, var5, var95, p50, p5, p95};
        
    }

    private static double[] getTrioSummaryStatistics(
            double childBeta,
            double motherBeta,
            double fatherBeta,
            double childSe,
            double motherSe,
            double fatherSe,
            double mafChild,
            double mafMother,
            double mafFather
    ) {

        double minMafMother = mafMother > 0.5 ? 1 - mafMother : mafMother;
        double minMafFather = mafFather > 0.5 ? 1 - mafFather : mafFather;
        
        ArrayList<Double> motherGenos = new ArrayList<>(nTriosSimulation);
        ArrayList<Double> fatherGenos = new ArrayList<>(nTriosSimulation);
        
        for (int i = 0 ; i < nTriosSimulation ; i++) {
            
            motherGenos.add(0.0);
            fatherGenos.add(0.0);
            
        }
        
        for (int i = 0 ; i < minMafMother * nTriosSimulation ; i++) {
            
            motherGenos.set(i, motherGenos.get(i) + 1);
            
        }
        
        Collections.shuffle(motherGenos);
        
        for (int i = 0 ; i < minMafMother * nTriosSimulation ; i++) {
            
            motherGenos.set(i, motherGenos.get(i) + 1);
            
        }
        
        Collections.shuffle(motherGenos);
        
        
        
        for (int i = 0 ; i < minMafFather * nTriosSimulation ; i++) {
            
            fatherGenos.set(i, fatherGenos.get(i) + 1);
            
        }
        
        Collections.shuffle(fatherGenos);
        
        for (int i = 0 ; i < minMafFather * nTriosSimulation ; i++) {
            
            fatherGenos.set(i, fatherGenos.get(i) + 1);
            
        }
        
        Collections.shuffle(fatherGenos);
        
        ArrayList<Double> phenos = new ArrayList<>(nTriosSimulation);
        double squaredError = 0.0;
        
        Random random = new Random();
        
        for (int i = 0 ; i < nTriosSimulation ; i++) {
            
            double motherGenotype = motherGenos.get(i);
            double fatherGenotype = fatherGenos.get(i);
            double childGenotype = (motherGenotype + fatherGenotype) / 2;
            
            double motherEffect = (motherBeta + motherSe * random.nextGaussian()) * motherGenotype;
            double fatherEffect = (fatherBeta + fatherSe * random.nextGaussian()) * fatherGenotype;
            double childEffect = (childBeta + childSe * random.nextGaussian()) * childGenotype;
            
            double phenoValue = random.nextGaussian() + motherEffect + fatherEffect + childEffect;
            
            phenos.add(phenoValue);
            
            double predictedPhenoValue = motherBeta * motherGenotype + fatherBeta * fatherGenotype + childEffect * childGenotype;
            
            squaredError += (predictedPhenoValue - phenoValue) * (predictedPhenoValue - phenoValue);
            
        }
        
        squaredError /= nTriosSimulation;
        
        double mean = phenos.stream().mapToDouble(a -> a).sum() / nTriosSimulation;
        
        double rss0 = 0.0;
        
        for (int i = 0 ; i < nTriosSimulation ; i++) {
            
            double phenoDist = phenos.get(i) - mean;
            rss0 += phenoDist * phenoDist;
            
        }
        
        rss0 /= nTriosSimulation;
        
        double varianceExplained = (rss0 - squaredError) / rss0;

        double p = getModelSignificance(rss0, 1, squaredError, 4, nTriosPValue);
        
        return new double[]{varianceExplained, p};
        
    }

    /**
     * Returns the significance of increasing the complexity of a simple model,
     * model1, with p1 parameters to a more complex model, model2, with p2
     * parameters. p1 < p2.
     *
     * @param model1RSS the residual sum of squares for model 1
     * @param p1 the number of parameters of model 1
     * @param model2RSS the residual sum of squares for model 2
     * @param p2 the number of parameters of model 2
     *
     * @return the significance
     */
    private static double getModelSignificance(
            double model1RSS,
            int p1,
            double model2RSS,
            int p2,
            int nValues
    ) {

        if (Double.isNaN(model1RSS) || Double.isNaN(model2RSS)) {

            return Double.NaN;

        }

        double numeratorDegreesOfFreedom = p2 - p1;
        double denominatorDegreesOfFreedom = nValues - p2;
        double x = ((model1RSS - model2RSS) / numeratorDegreesOfFreedom) / (model2RSS / denominatorDegreesOfFreedom);

        FDistribution fDistribution = new FDistribution(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom);

        return 1.0 - fDistribution.cumulativeProbability(x);

    }

}
