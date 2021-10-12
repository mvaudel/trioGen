package no.uib.triogen.scripts_marc.misc;

import java.io.File;
import java.util.HashSet;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 *
 *
 * @author Marc Vaudel
 */
public class GetTargetHits {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        HashSet<String> h2020Targets = getIds(new File("C:\\Github\\SFF21\\resources\\targets\\targets_h2020"));
        HashSet<String> pwTargets = getIds(new File("C:\\Github\\SFF21\\resources\\targets\\targets_pw"));

        try (SimpleFileWriter writer = new SimpleFileWriter(new File("C:\\Github\\SFF21\\resources\\top_hits"), false)) {

            writer.writeLine("phenotype", "snp", "chr", "pos", "tested_allele", "other_allele", "beta", "se", "p");

            // PW
            System.out.println("Processing PW");

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(new File("C:\\Projects\\EGG\\placenta_weight\\meta\\pw_fetal_sex_gest_clean.gz"))) {

                String line = reader.readLine();

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String rsid = lineSplit[0];

                    if (h2020Targets.contains(rsid) || pwTargets.contains(rsid)) {

                        String chr = lineSplit[1];
                        String pos = lineSplit[2];
                        String tested_allele = lineSplit[3];
                        String other_allele = lineSplit[4];
                        String beta = lineSplit[9];
                        String se = lineSplit[10];
                        String p = lineSplit[11];

                        writer.writeLine("PW", rsid, chr, pos, tested_allele, other_allele, beta, se, p);

                    }
                }
            }

            // H2020
            String[] filePatterns = new String[]{
                "C:\\Projects\\h2020\\gwas\\z_bmiageI-stats-bgen.gz",
                "C:\\Projects\\h2020\\gwas\\z_bmiageI-stats-bgen-chrX.gz"
            };

            for (int ageI = 0; ageI <= 11; ageI++) {

                System.out.println("Processing H2020 " + ageI);

                for (String filePattern : filePatterns) {

                    String filePath = filePattern.replaceAll("ageI", Integer.toString(ageI));

                    try (SimpleFileReader reader = SimpleFileReader.getFileReader(new File(filePath))) {

                        String line = reader.readLine();

                        while ((line = reader.readLine()) != null) {

                            String[] lineSplit = line.split("\t");

                            String rsid = lineSplit[0];

                            if (h2020Targets.contains(rsid) || pwTargets.contains(rsid)) {

                                String chr = lineSplit[1];
                                String pos = lineSplit[2];
                                String tested_allele = lineSplit[4];
                                String other_allele = lineSplit[5];
                                String beta = lineSplit[10];
                                String se = lineSplit[11];
                                String p = lineSplit[15];

                                writer.writeLine(Integer.toString(ageI), rsid, chr, pos, tested_allele, other_allele, beta, se, p);

                            }
                        }
                    }
                }
            }
        }
    }

    public static HashSet<String> getIds(File file) {

        HashSet<String> result = new HashSet<>();

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(file)) {

            String line = reader.readLine();

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split("\t");

                result.add(lineSplit[0]);

            }
        }

        return result;

    }

}
