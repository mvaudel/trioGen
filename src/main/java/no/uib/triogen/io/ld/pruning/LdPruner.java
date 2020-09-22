package no.uib.triogen.io.ld.pruning;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.utils.Utils;

/**
 * Simple class for LD pruning.
 *
 * @author Marc Vaudel
 */
public class LdPruner {

    /**
     * The ld matrix file path.
     */
    private final String ldMatrixFilePath;
    /**
     * The results file to prune.
     */
    private final File resultsFile;
    /**
     * File where to write the output.
     */
    private final File destinationFile;
    /**
     * The minimal r2 to export.
     */
    private final double minR2;
    /**
     * The maximal p to export.
     */
    private final double maxP;
    /**
     * The name of the variant id column.
     */
    private final String variantIdColName;
    /**
     * The name of the p-value column.
     */
    private final String pColName;
    /**
     * The name of the pheno column.
     */
    private final String phenoColName;
    /**
     * The name of the contig column.
     */
    private final String contigColName;
    /**
     * The column separator.
     */
    private final String separator;

    private LdMatrixReader defaultLdMatrixReader;

    public LdPruner(
            String ldMatrixFilePath,
            File resultsFile,
            File destinationFile,
            double minR2,
            double maxP,
            String pColName,
            String variantIdColName,
            String phenoColName,
            String contigColName,
            String separator
    ) {

        this.ldMatrixFilePath = ldMatrixFilePath;
        this.resultsFile = resultsFile;
        this.destinationFile = destinationFile;
        this.minR2 = minR2;
        this.maxP = maxP;
        this.variantIdColName = variantIdColName;
        this.pColName = pColName;
        this.phenoColName = phenoColName;
        this.contigColName = contigColName;
        this.separator = separator;

    }

    public void prune() {

        Instant begin = Instant.now();

        System.out.println("Getting p-values from " + resultsFile.getAbsolutePath() + ".");

        HashMap<String, HashMap<String, TreeMap<Double, TreeMap<String, String>>>> contigPhenoPvalueVariantLineMap = new HashMap<>();

        int lineNumber = 0;

        try ( SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            int variantIdColIndex = -1;
            int pColIndex = -1;
            int phenoColIndex = -1;
            int contigColIndex = -1;

            String headerLine = reader.readLine();

            lineNumber++;

            String[] lineSplit = headerLine.split(separator);

            for (int i = 0; i < lineSplit.length; i++) {

                if (lineSplit[i].equals(variantIdColName)) {

                    variantIdColIndex = i;

                }
                if (lineSplit[i].equals(pColName)) {

                    pColIndex = i;

                }
                if (phenoColName != null && lineSplit[i].equals(phenoColName)) {

                    phenoColIndex = i;

                }
                if (lineSplit[i].equals(contigColName)) {

                    contigColIndex = i;

                }
            }

            if (variantIdColIndex == -1) {

                throw new IllegalArgumentException("Variant id column '" + variantIdColName + "' not found in " + resultsFile + ".\n" + headerLine);

            }
            if (pColIndex == -1) {

                throw new IllegalArgumentException("P-value column '" + pColName + "' not found in " + resultsFile + ".\n" + headerLine);

            }
            if (phenoColName != null && phenoColIndex == -1) {

                throw new IllegalArgumentException("Pheno column '" + phenoColName + "' not found in " + resultsFile + ".\n" + headerLine);

            }
            if (contigColIndex == -1) {

                throw new IllegalArgumentException("contig column '" + contigColName + "' not found in " + resultsFile + ".\n" + headerLine);

            }

            String line;

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(separator);

                String contig = lineSplit[contigColIndex];
                String pheno = phenoColIndex == -1 ? "dummy" : lineSplit[phenoColIndex];
                String variantId = lineSplit[variantIdColIndex];
                String pValueString = lineSplit[pColIndex];

                double pValue;

                try {

                    pValue = Double.parseDouble(pValueString);

                } catch (Exception e) {

                    throw new IllegalArgumentException("P-value '" + pValueString + "' at line " + lineNumber + " in " + resultsFile + " could not be parsed as a number.");

                }

                if (pValue <= maxP) {

                    HashMap<String, TreeMap<Double, TreeMap<String, String>>> phenoPvalueVariantLineMap = contigPhenoPvalueVariantLineMap.get(pheno);

                    if (phenoPvalueVariantLineMap == null) {

                        phenoPvalueVariantLineMap = new HashMap<>();
                        contigPhenoPvalueVariantLineMap.put(pheno, phenoPvalueVariantLineMap);

                    }

                    TreeMap<Double, TreeMap<String, String>> pvalueVariantLineMap = phenoPvalueVariantLineMap.get(contig);

                    if (pvalueVariantLineMap == null) {

                        pvalueVariantLineMap = new TreeMap<>();
                        phenoPvalueVariantLineMap.put(contig, pvalueVariantLineMap);

                    }

                    TreeMap<String, String> variantLineMap = pvalueVariantLineMap.get(pValue);

                    if (variantLineMap == null) {

                        variantLineMap = new TreeMap<>();
                        pvalueVariantLineMap.put(pValue, variantLineMap);

                    }

                    variantLineMap.put(variantId, line);

                }

                lineNumber++;

            }

            Instant end = Instant.now();

            long timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("Getting p-values finished (" + timeInSec + " s)");

            begin = Instant.now();

            System.out.println("LD pruning.");

            TreeMap<String, TreeMap<String, ArrayList<String>>> prunedLines = prune(contigPhenoPvalueVariantLineMap);

            end = Instant.now();

            timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("Pruning finished (" + timeInSec + " s)");

            begin = Instant.now();

            System.out.println("Exporting to " + destinationFile + ".");

            try ( SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

                writer.writeLine(headerLine);

                prunedLines.values().stream()
                        .flatMap(
                                submap -> submap.values().stream()
                        )
                        .flatMap(
                                lines -> lines.stream()
                        )
                        .forEach(
                                prunedLine -> writer.writeLine(prunedLine)
                        );
            }

            end = Instant.now();

            timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("Export finished (" + timeInSec + " s)");

        }
    }

    private TreeMap<String, TreeMap<String, ArrayList<String>>> prune(
            HashMap<String, HashMap<String, TreeMap<Double, TreeMap<String, String>>>> contigMap
    ) {

        return contigMap.entrySet()
                .parallelStream()
                .collect(
                        Collectors.toMap(
                                Entry::getKey,
                                entry -> prune(entry.getValue(), getLdMatrixReader(entry.getKey())),
                                (a, b) -> a,
                                TreeMap::new
                        )
                );

    }

    private LdMatrixReader getLdMatrixReader(
            String contig
    ) {

        if (contigColName != null) {

            String contigLdMatrixFilePath = ldMatrixFilePath.replace(Utils.CONTIG_WILDCARD, contig);

            try {

                return new LdMatrixReader(new File(contigLdMatrixFilePath));

            } catch (Exception e) {

                throw new RuntimeException(e);

            }

        } else {

            return defaultLdMatrixReader;

        }
    }

    private TreeMap<String, ArrayList<String>> prune(
            HashMap<String, TreeMap<Double, TreeMap<String, String>>> phenoMap,
            LdMatrixReader ldMatrixReader
    ) {

        return phenoMap.entrySet()
                .parallelStream()
                .collect(
                        Collectors.toMap(
                                Entry::getKey,
                                entry -> prune(entry.getValue(), ldMatrixReader),
                                (a, b) -> a,
                                TreeMap::new
                        )
                );

    }

    private ArrayList<String> prune(
            TreeMap<Double, TreeMap<String, String>> pValuesMap,
            LdMatrixReader ldMatrixReader
    ) {

        try {

            ArrayList<String> result = new ArrayList<>();

            HashSet<String> inspectedSnp = new HashSet<>();

            for (Entry<Double, TreeMap<String, String>> entry1 : pValuesMap.entrySet()) {

                for (Entry<String, String> entry2 : entry1.getValue().entrySet()) {

                    String variantId = entry2.getKey();

                    if (!inspectedSnp.contains(variantId)) {

                        result.add(entry2.getValue());
                        inspectedSnp.add(variantId);

                    }

                    HashMap<String, Double> variantLdMap = ldMatrixReader.getR2(variantId);

                    if (variantLdMap != null) {

                        for (Entry<String, Double> entry3 : variantLdMap.entrySet()) {

                            if (entry3.getValue() >= minR2) {

                                inspectedSnp.add(entry3.getKey());

                            }
                        }
                    }
                }
            }

            return result;

        } catch (Exception e) {

            throw new RuntimeException(e);

        }
    }
}
