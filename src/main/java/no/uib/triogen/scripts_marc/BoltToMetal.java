package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * Converts BOLT-LMM results to metal
 *
 * @author Marc Vaudel
 */
public class BoltToMetal {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String boltResultsFolder = "/mnt/work/marc/moba/run/bolt/bolt_output";
        String metaFolder = "/mnt/work/marc/moba/pwbw/prs/meta/moba_training";

        String[] phenotypes = new String[]{"z_placenta_weight", "z_weight0", "z_pw_bw_ratio"};
        String[] prefixes = new String[]{"child", "mother", "father"};

        for (String phenotype : phenotypes) {
            for (String prefix : prefixes) {

                System.out.println(Instant.now() + "    Processing GWAS on " + phenotype + " against " + prefix + " genome.");

                Instant begin = Instant.now();

                File metaFile = new File(metaFolder + "/" + prefix + "Geno_" + phenotype + "_meta.gz");

                try (SimpleFileWriter writer = new SimpleFileWriter(metaFile, true)) {

                    writer.writeLine(
                            "SNP", "SNPID", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "N", "EAF", "BETA", "SE", "PVAL"
                    );

                    File boltResultsFile = new File(boltResultsFolder + "/" + phenotype + "/" + prefix + "Geno_" + phenotype + "-stats-bgen_training.gz");

                    if (!boltResultsFile.exists()) {

                        throw new IllegalArgumentException("Bolt result file not found " + boltResultsFile + ".");

                    }

                    File logFile = new File("/mnt/work/marc/moba/run/bolt/bolt_output/" + phenotype + "/" + prefix + "Geno_" + phenotype + "-runlog_training.log");

                    if (!logFile.exists()) {

                        throw new IllegalArgumentException("Bolt log file not found " + logFile + ".");

                    }

                    int nSamples = getNSamples(logFile);

                    process(boltResultsFile, nSamples, writer);

                    boltResultsFile = new File(boltResultsFolder + "/" + phenotype + "/" + prefix + "Geno_" + phenotype + "-stats-bgen-chrX_training.gz");

                    if (!boltResultsFile.exists()) {

                        throw new IllegalArgumentException("Bolt result file for chr X not found " + boltResultsFile + ".");

                    }

                    logFile = new File("/mnt/work/marc/moba/run/bolt/bolt_output/" + phenotype + "/" + prefix + "Geno_" + phenotype + "-runlog-chrX_training.log");

                    if (!logFile.exists()) {

                        throw new IllegalArgumentException("Bolt log file not found " + logFile + ".");

                    }

                    nSamples = getNSamples(logFile);

                    process(boltResultsFile, nSamples, writer);

                    Instant end = Instant.now();

                    long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                    System.out.println(Instant.now() + "    Processing GWAS on " + phenotype + " against " + prefix + " genome finished (" + durationSeconds + " s)");

                }
            }
        }
    }

    private final static String nUsed = "Number of individuals used in analysis: Nused = ";

    private static int getNSamples(File logFile) {

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(logFile, false)) {

            String line;
            while ((line = reader.readLine()) != null) {

                if (line.startsWith(nUsed)) {

                    String stringValue = line.substring(nUsed.length()).trim();

                    return Integer.parseInt(stringValue);

                }
            }
        }

        throw new IllegalArgumentException("Number of individuals not found in log file " + logFile + ".");

    }

    private final static String headerLine = "SNP\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tINFO\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\tCHISQ_BOLT_LMM_INF\tP_BOLT_LMM_INF\tCHISQ_BOLT_LMM\tP_BOLT_LMM";

    private static void process(
            File boltResultsFile,
            int nSamples,
            SimpleFileWriter writer
    ) {

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(boltResultsFile, false)) {

            String line = reader.readLine();

            if (!line.equals(headerLine)) {

                throw new IllegalArgumentException("Unexpected header in " + boltResultsFile + ".\n" + line);

            }

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                String rsid = lineSplit[0];
                String chr = lineSplit[1];
                String bp = lineSplit[2];
                String testedAllele = lineSplit[4];
                String otherAllele = lineSplit[5];
                String testedAlleleFreq = lineSplit[6];
                String beta = lineSplit[10];
                String se = lineSplit[11];
                String p = lineSplit[15];

                TreeSet<String> alleles = new TreeSet<>();
                alleles.add(testedAllele);
                alleles.add(otherAllele);

                String snp = chr + ":" + bp + "_" + alleles.stream().collect(Collectors.joining("_"));

                writer.writeLine(
                        snp, rsid, testedAllele, otherAllele, Integer.toString(nSamples), testedAlleleFreq, beta, se, p
                );

            }
        }
    }
}
