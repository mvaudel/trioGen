/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uib.triogen.scripts_marc;

import java.io.File;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.cell_rk.utils.SimpleFileWriter;

/**
 *
 * @author mvaudel
 */
public class DebugPw {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        HashMap<String, Double> huntFatherBeta = new HashMap<>();

        File resultsFile = new File("C:\\Github\\placenta_weight\\tmp\\EGG_HRC_BW6.PW.father.sex_gest.HUNT.european.CF.20200702.txt.gz");
        File targetsFile = new File("C:\\Github\\placenta_weight\\tmp\\father_hunt_targets");

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(targetsFile, false)) {

                String line = reader.readLine();

                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");
                    

                    String contig = lineSplit[4];
                    String pos = lineSplit[5];
                    String alleleA = lineSplit[6];
                    String alleleB = lineSplit[7];
                    double maf = Double.parseDouble(lineSplit[12]);
                    double beta = Double.parseDouble(lineSplit[15]);

                    if (alleleA.compareTo(alleleB) > 0) {

                        alleleB = alleleA;
                        beta = -beta;

                    }

                    String key = String.join("\t", contig, pos, alleleA, alleleB);
                    
                    if (maf >= 0.1 && maf <= 0.9) {

                    huntFatherBeta.put(key, beta);
                    
                    }

                    if (contig.equals("1") && pos.equals("155001281")) {

                        writer.writeLine(line);

                    }
                    if (contig.equals("5") && pos.equals("132444128")) {

                        writer.writeLine(line);

                    }
                }
            }
        }

        HashMap<String, Double> mobaFatherBeta = new HashMap<>(huntFatherBeta.size());

        resultsFile = new File("C:\\Github\\placenta_weight\\tmp\\EGG_HRC_BW6.PW.father.sex_gest.MoBa.european.CF.20200623.txt.gz");
        targetsFile = new File("C:\\Github\\placenta_weight\\tmp\\father_moba_targets");

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(targetsFile, false)) {

                String line = reader.readLine();

                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String contig = lineSplit[4];
                    String pos = lineSplit[5];
                    String alleleA = lineSplit[6];
                    String alleleB = lineSplit[7];
                    double maf = Double.parseDouble(lineSplit[12]);
                    double beta = Double.parseDouble(lineSplit[15]);

                    if (alleleA.compareTo(alleleB) > 0) {

                        alleleB = alleleA;
                        beta = -beta;

                    }

                    String key = String.join("\t", contig, pos, alleleA, alleleB);

                    if (huntFatherBeta.containsKey(key) && maf >= 0.1 && maf <= 0.9) {

                        mobaFatherBeta.put(key, beta);

                    }

                    if (contig.equals("1") && pos.equals("155001281")) {

                        writer.writeLine(line);

                    }
                    if (contig.equals("5") && pos.equals("132444128")) {

                        writer.writeLine(line);

                    }
                }
            }
        }

        HashMap<String, Double> huntMotherBeta = new HashMap<>(huntFatherBeta.size());

        resultsFile = new File("C:\\Github\\placenta_weight\\tmp\\EGG_HRC_BW6.PW.mother.sex_gest.HUNT.european.CF.20200702.txt.gz");
        targetsFile = new File("C:\\Github\\placenta_weight\\tmp\\mother_hunt_targets");

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(targetsFile, false)) {

                String line = reader.readLine();

                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String contig = lineSplit[4];
                    String pos = lineSplit[5];
                    String alleleA = lineSplit[6];
                    String alleleB = lineSplit[7];
                    double maf = Double.parseDouble(lineSplit[12]);
                    double beta = Double.parseDouble(lineSplit[15]);

                    if (alleleA.compareTo(alleleB) > 0) {

                        alleleB = alleleA;
                        beta = -beta;

                    }

                    String key = String.join("\t", contig, pos, alleleA, alleleB);

                    if (mobaFatherBeta.containsKey(key) && maf >= 0.1 && maf <= 0.9) {

                        huntMotherBeta.put(key, beta);

                    }

                    if (contig.equals("1") && pos.equals("155001281")) {

                        writer.writeLine(line);

                    }
                    if (contig.equals("5") && pos.equals("132444128")) {

                        writer.writeLine(line);

                    }
                }
            }
        }

        HashMap<String, Double> mobaMotherBeta = new HashMap<>(huntMotherBeta.size());

        resultsFile = new File("C:\\Github\\placenta_weight\\tmp\\EGG_HRC_BW6.PW.mother.sex_gest.MoBa.european.CF.20200623.txt.gz");
        targetsFile = new File("C:\\Github\\placenta_weight\\tmp\\mother_moba_targets");

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(targetsFile, false)) {

                String line = reader.readLine();

                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String contig = lineSplit[4];
                    String pos = lineSplit[5];
                    String alleleA = lineSplit[6];
                    String alleleB = lineSplit[7];
                    double maf = Double.parseDouble(lineSplit[12]);
                    double beta = Double.parseDouble(lineSplit[15]);

                    if (alleleA.compareTo(alleleB) > 0) {

                        alleleB = alleleA;
                        beta = -beta;

                    }

                    String key = String.join("\t", contig, pos, alleleA, alleleB);

                    if (huntMotherBeta.containsKey(key) && maf >= 0.1 && maf <= 0.9) {

                        mobaMotherBeta.put(key, beta);

                    }

                    if (contig.equals("1") && pos.equals("155001281")) {

                        writer.writeLine(line);

                    }
                    if (contig.equals("5") && pos.equals("132444128")) {

                        writer.writeLine(line);

                    }
                }
            }
        }

        try (SimpleFileWriter writer = new SimpleFileWriter(new File("C:\\Github\\placenta_weight\\tmp\\beta_father"), true)) {

            writer.writeLine(
                    String.join("\t",
                            "chr",
                            "bp",
                            "a",
                            "b",
                            "beta_mother_hunt",
                            "beta_mother_moba",
                            "beta_father_hunt",
                            "beta_father_moba"
                    )
            );

            for (Entry<String, Double> entry : mobaMotherBeta.entrySet()) {

                String key = entry.getKey();
                Double betaMobaMother = entry.getValue();
                Double betaHuntMother = huntMotherBeta.get(key);
                Double betaMobaFather = mobaFatherBeta.get(key);
                Double betaHuntFather = huntFatherBeta.get(key);

                writer.writeLine(
                        String.join("\t",
                                key,
                                Double.toString(betaHuntMother),
                                Double.toString(betaMobaMother),
                                Double.toString(betaHuntFather),
                                Double.toString(betaMobaFather)
                        )
                );
            }
        }
    }

    private static String swapBeta(String beta) {

        if (beta.charAt(0) == '-') {

            return beta.substring(1);

        } else {

            StringBuilder sb = new StringBuilder(beta.length() + 1);

            return sb.append('-')
                    .append(beta)
                    .toString();

        }

    }

}
