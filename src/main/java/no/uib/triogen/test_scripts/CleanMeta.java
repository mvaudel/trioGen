package no.uib.triogen.test_scripts;

import java.io.File;
import java.time.Instant;
import java.util.HashMap;
import java.util.TreeSet;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Convenience script cleaning the results of the pw meta.
 *
 * @author Marc Vaudel
 */
public class CleanMeta {

    public static final double MAC_THRESHOLD = 10;
    public static final double INFO_THRESHOLD = 0.4;
    public static final double ROUNDING = 1e-6;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            File folder = new File(args[0]);

            for (File file : folder.listFiles()) {

                if (file.getName().endsWith(".gz")) {

                    processFile(file);

                }

            }
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

    private static void processFile(File file) {

        String filePath = file.getAbsolutePath();

        Instant start = Instant.now();

        System.out.println(start + "    Processing " + filePath);

        File outFile = new File(new File(file.getParentFile(), "clean"), file.getName());

        try ( SimpleFileReader reader = SimpleFileReader.getFileReader(file)) {

            try ( SimpleFileWriter writer = new SimpleFileWriter(outFile, true)) {

                String line = reader.readLine();

                String separator = "\t";

                String[] lineSplit = line.split(separator);

                if (lineSplit.length == 1) {

                    separator = " ";
                    lineSplit = line.split(separator);

                    if (lineSplit.length == 1) {

                        throw new IllegalArgumentException("Separator not found for file " + file + ".");

                    }

                    line = line.replace(' ', '\t');

                }

                line = String.join("\t", "SNP", line);

                writer.writeLine(line);

                HashMap<String, Integer> columnMap = new HashMap<>(lineSplit.length);

                for (int i = 0; i < lineSplit.length; i++) {

                    columnMap.put(lineSplit[i], i);

                }

                int lineNumber = 1;

                while ((line = reader.readLine()) != null) {

                    // Read new line, sanity check on the number of clumn
                    lineNumber++;

                    lineSplit = line.split(separator);

                    if (lineSplit.length != columnMap.size()) {

                        throw new IllegalArgumentException(lineSplit.length + " columns found at line " + lineNumber + " in file " + file + " where " + columnMap.size() + " expected.");

                    }

                    // Read info, skip low values
                    double info;

                    if (columnMap.containsKey("INFO")) {

                        String cellValue = lineSplit[columnMap.get("INFO")];

                        if (cellValue.equals(".") || cellValue.equals("NA") || cellValue.equals("NaN")) {

                            info = 1.0;

                        } else {

                            try {

                                info = Double.parseDouble(cellValue);

                            } catch (Exception e) {

                                throw new IllegalArgumentException("info " + cellValue + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                            }

                            if (info < 0 && info >= -ROUNDING) {

                                info = 0;

                            }
                            if (info > 1.0 && info <= 1.0 + ROUNDING) {

                                info = 1.0;

                            }

                            if (info < 0 || info > 1.0) {

                                throw new IllegalArgumentException("info " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                            }
                        }
                    } else if (columnMap.containsKey("INFO_fem") || columnMap.containsKey("INFO_male")) {

                        double infoFemales;

                        String cellValue = lineSplit[columnMap.get("INFO_fem")];

                        if (cellValue.equals(".") || cellValue.equals("NA") || cellValue.equals("NaN")) {

                            infoFemales = 1.0;

                        } else {

                            try {

                                infoFemales = Double.parseDouble(cellValue);

                            } catch (Exception e) {

                                throw new IllegalArgumentException("infoFemales " + cellValue + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                            }

                            if (infoFemales < 0 && infoFemales >= -ROUNDING) {

                                infoFemales = 0;

                            }
                            if (infoFemales > 1.0 && infoFemales <= 1.0 + ROUNDING) {

                                infoFemales = 1.0;

                            }

                            if (infoFemales < 0 || infoFemales > 1.0) {

                                throw new IllegalArgumentException("infoFemales " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                            }
                        }

                        cellValue = lineSplit[columnMap.get("INFO_male")];

                        double infoMales;

                        if (cellValue.equals(".") || cellValue.equals("NA") || cellValue.equals("NaN")) {

                            infoMales = 1.0;

                        } else {

                            try {

                                infoMales = Double.parseDouble(cellValue);

                            } catch (Exception e) {

                                throw new IllegalArgumentException("infoMales " + cellValue + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                            }

                            if (infoMales < 0 && infoMales >= -ROUNDING) {

                                infoMales = 0;

                            }
                            if (infoMales > 1.0 && infoMales <= 1.0 + ROUNDING) {

                                infoMales = 1.0;

                            }

                            if (infoMales < 0 || infoMales > 1.0) {

                                throw new IllegalArgumentException("infoMales " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                            }
                        }

                        info = (infoFemales + infoMales) / 2;

                    } else {

                        throw new IllegalArgumentException("info column not found in file " + file + ".");

                    }

                    if (info < INFO_THRESHOLD) {

                        continue;

                    }

                    // Read se, skip NAs
                    String cellValue = lineSplit[columnMap.get("SE")];

                    if (cellValue.equals("NA") || cellValue.equals("NaN") || cellValue.equals(".")) {

                        continue;

                    }

                    double se;

                    try {

                        se = Double.parseDouble(cellValue);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("se " + cellValue + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                    }

                    if (se <= 0) {

                        throw new IllegalArgumentException("se " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                    }

                    // Read maf, skip low allele count
                    cellValue = lineSplit[columnMap.get("EAF")];
                    double maf;

                    try {

                        maf = Double.parseDouble(cellValue);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("maf " + cellValue + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                    }

                    if (maf < 0 || maf > 1.0) {

                        throw new IllegalArgumentException("maf " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                    }

                    if (maf > 0.5) {

                        maf = 1.0 - maf;

                    }

                    cellValue = lineSplit[columnMap.get("N")];
                    double n;

                    try {

                        n = Integer.parseInt(cellValue);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("N " + cellValue + " could not be parsed as an integer at line " + lineNumber + " in file " + file + ".");

                    }

                    if (n <= 0) {

                        throw new IllegalArgumentException("N " + cellValue + " out of range at line " + lineNumber + " in file " + file + ".");

                    }

                    double mac = 2.0 * n * maf;

                    if (mac < MAC_THRESHOLD) {

                        continue;

                    }

                    // Parse contig, use 23 for the 'X' chromosome
                    String contig = lineSplit[columnMap.get("CHR")];

                    if (contig.equals("0X") || contig.equals("X") || contig.equals("NA")) {

                        contig = "23";

                    }

                    try {

                        Double.parseDouble(contig);

                    } catch (Exception e) {

                        throw new IllegalArgumentException("Contig " + contig + " could not be parsed as a number at line " + lineNumber + " in file " + file + ".");

                    }

                    // Set snp id as CHR:POS_A_B where A and B are in alphabetical order
                    TreeSet<String> alleles = new TreeSet<>();

                    alleles.add(lineSplit[columnMap.get("EFFECT_ALLELE")]);
                    alleles.add(lineSplit[columnMap.get("NON_EFFECT_ALLELE")]);

                    String bp = lineSplit[columnMap.get("POS")];

                    StringBuilder snpIdBuilder = new StringBuilder()
                            .append(contig)
                            .append(':')
                            .append(bp);

                    alleles.forEach(
                            allele -> snpIdBuilder
                                    .append('_')
                                    .append(allele)
                    );

                    String snpId = snpIdBuilder.toString();

                    // Make sure the line is tab separated and add the snp id in first column
                    if (separator.equals(" ")) {

                        line = line.replace(' ', '\t');

                    }

                    line = String.join("\t", snpId, line);

                    // Write to the result file
                    writer.writeLine(line);

                }
            }
        }

        Instant end = Instant.now();

        long timeInSeconds = end.getEpochSecond() - start.getEpochSecond();

        System.out.println(end + "    Processing " + filePath + " (" + timeInSeconds + " s)");

    }

}
