package no.uib.triogen.processing.ld.pruning;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.model.annotation.ProxyCoordinates;
import no.uib.triogen.model.annotation.ensembl.EnsemblAPI;
import no.uib.triogen.model.annotation.ld_link.LDproxy;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.utils.Utils;

/**
 * Simple LD pruning.
 *
 * @author Marc Vaudel
 */
public class SimpleLdPruner {

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
     * The name of the rsid column.
     */
    private final String rsidColName;
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
    /**
     * The LD matrix file reader.
     */
    private LdMatrixReader defaultLdMatrixReader;
    /**
     * The build number to use when querying Ensembl.
     */
    private final int buildNumber;
    /**
     * The reference population to use for Ensembl.
     */
    private final String ensemblPopulation;
    /**
     * The reference population to use for LDlink.
     */
    private final String ldLinkPopulation;
    /**
     * The token to use for LDlink.
     */
    private final String ldLinkToken;

    /**
     * Constructor.
     *
     * @param ldMatrixFilePath The path to the LD matrices.
     * @param resultsFile The results file to prune.
     * @param destinationFile The file where to write the results.
     * @param minR2 The minimal R2 to consider that two hits cannot be
     * considered independent.
     * @param maxP The maximal p-value to consider.
     * @param pColName The p-value column name.
     * @param variantIdColName The variant id column name.
     * @param rsidColName The rsid column name.
     * @param phenoColName The phenotype column name.
     * @param contigColName The contig column name.
     * @param separator The separator for the columns.
     * @param buildNumber The build number to use when querying Ensembl.
     * @param ensemblPopulation The reference population to use for Ensembl.
     * @param ldLinkPopulation The reference population to use for LDlink.
     * @param ldLinkToken The token to use for LDlink.
     */
    public SimpleLdPruner(
            String ldMatrixFilePath,
            File resultsFile,
            File destinationFile,
            double minR2,
            double maxP,
            String pColName,
            String variantIdColName,
            String rsidColName,
            String phenoColName,
            String contigColName,
            String separator,
            int buildNumber,
            String ensemblPopulation,
            String ldLinkPopulation,
            String ldLinkToken
    ) {

        this.ldMatrixFilePath = ldMatrixFilePath;
        this.resultsFile = resultsFile;
        this.destinationFile = destinationFile;
        this.minR2 = minR2;
        this.maxP = maxP;
        this.variantIdColName = variantIdColName;
        this.rsidColName = rsidColName;
        this.pColName = pColName;
        this.phenoColName = phenoColName;
        this.contigColName = contigColName;
        this.separator = separator;
        this.buildNumber = buildNumber;
        this.ensemblPopulation = ensemblPopulation;
        this.ldLinkPopulation = ldLinkPopulation;
        this.ldLinkToken = ldLinkToken;

    }

    /**
     * Runs the pruning.
     */
    public void run() {

        Instant begin = Instant.now();

        System.out.println("Getting p-values from " + resultsFile.getAbsolutePath() + ".");

        HashMap<String, HashMap<String, TreeMap<Double, TreeMap<String, String>>>> contigPhenoPvalueVariantLineMap = new HashMap<>();
        HashMap<String, String> variantToRsidMap = new HashMap<>();

        int lineNumber = 0;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultsFile)) {

            int variantIdColIndex = -1;
            int rsidColIndex = -1;
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
                if (lineSplit[i].equals(rsidColName)) {

                    rsidColIndex = i;

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
            if ((ensemblPopulation != null || ldLinkPopulation != null) && rsidColIndex == -1) {

                throw new IllegalArgumentException("Rsid column '" + variantIdColName + "' not found in " + resultsFile + ".\n" + headerLine);

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
                String rsid = rsidColIndex == -1 ? variantId : lineSplit[rsidColIndex];
                String pValueString = lineSplit[pColIndex];

                double pValue;

                try {

                    pValue = Double.parseDouble(pValueString);

                } catch (Exception e) {

                    throw new IllegalArgumentException("P-value '" + pValueString + "' at line " + lineNumber + " in " + resultsFile + " could not be parsed as a number.");

                }

                if (pValue <= maxP) {

                    variantToRsidMap.put(variantId, rsid);

                    HashMap<String, TreeMap<Double, TreeMap<String, String>>> phenoPvalueVariantLineMap = contigPhenoPvalueVariantLineMap.get(contig);

                    if (phenoPvalueVariantLineMap == null) {

                        phenoPvalueVariantLineMap = new HashMap<>();
                        contigPhenoPvalueVariantLineMap.put(contig, phenoPvalueVariantLineMap);

                    }

                    TreeMap<Double, TreeMap<String, String>> pvalueVariantLineMap = phenoPvalueVariantLineMap.get(pheno);

                    if (pvalueVariantLineMap == null) {

                        pvalueVariantLineMap = new TreeMap<>();
                        phenoPvalueVariantLineMap.put(pheno, pvalueVariantLineMap);

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

            TreeMap<String, TreeMap<String, ArrayList<String>>> prunedLines = prune(contigPhenoPvalueVariantLineMap, variantToRsidMap);

            end = Instant.now();

            timeInSec = end.getEpochSecond() - begin.getEpochSecond();

            System.out.println("Pruning finished (" + timeInSec + " s)");

            begin = Instant.now();

            System.out.println("Exporting to " + destinationFile + ".");

            try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

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

    /**
     * Prunes a map of p-values.
     *
     * @param contigMap The map of p-values per contig.
     * @param variantToRsidMap The variant id to rsid map.
     *
     * @return The pruned map.
     */
    private TreeMap<String, TreeMap<String, ArrayList<String>>> prune(
            HashMap<String, HashMap<String, TreeMap<Double, TreeMap<String, String>>>> contigMap,
            HashMap<String, String> variantToRsidMap
    ) {

        return contigMap.entrySet()
                .parallelStream()
                .filter(
                        entry -> getLdMatrixReader(entry.getKey()) != null
                )
                .collect(
                        Collectors.toMap(
                                Entry::getKey,
                                entry -> prune(entry.getValue(), variantToRsidMap, getLdMatrixReader(entry.getKey())),
                                (a, b) -> a,
                                TreeMap::new
                        )
                );

    }

    /**
     * Returns the LD matrix reader for the given contig.
     *
     * @param contig The contig name.
     *
     * @return The LD matrix reader for the given contig.
     */
    private LdMatrixReader getLdMatrixReader(
            String contig
    ) {

        if (contigColName != null) {

            String contigLdMatrixFilePath = ldMatrixFilePath.replace(Utils.CHROMOSOME_WILDCARD, contig);

            File contigLdFile = new File(contigLdMatrixFilePath);

            if (!contigLdFile.exists()) {

                return null;

            }

            try {

                return new LdMatrixReader(contigLdFile);

            } catch (Exception e) {

                throw new RuntimeException(e);

            }

        } else {

            return defaultLdMatrixReader;

        }
    }

    /**
     * Prunes the p-values for a phenotype.
     *
     * @param phenoMap The map of p-values for a given phenotype.
     * @param variantToRsidMap The variant id to rsid map.
     * @param ldMatrixReader The LD matrix reader.
     *
     * @return The pruned values for this phenotype.
     */
    private TreeMap<String, ArrayList<String>> prune(
            HashMap<String, TreeMap<Double, TreeMap<String, String>>> phenoMap,
            HashMap<String, String> variantToRsidMap,
            LdMatrixReader ldMatrixReader
    ) {

        return phenoMap.entrySet()
                .parallelStream()
                .collect(
                        Collectors.toMap(
                                Entry::getKey,
                                entry -> prune(entry.getValue(), variantToRsidMap, ldMatrixReader),
                                (a, b) -> a,
                                TreeMap::new
                        )
                );

    }

    /**
     * Prunes a p-value map.
     *
     * @param pValuesMap The p-values.
     * @param variantToRsidMap The variant id to rsid map.
     * @param ldMatrixReader The LD matrix reader.
     *
     * @return The pruned p-values.
     */
    private ArrayList<String> prune(
            TreeMap<Double, TreeMap<String, String>> pValuesMap,
            HashMap<String, String> variantToRsidMap,
            LdMatrixReader ldMatrixReader
    ) {

        try {

            ArrayList<String> result = new ArrayList<>();

            HashSet<String> inspectedSnp = new HashSet<>();

            for (Entry<Double, TreeMap<String, String>> entry1 : pValuesMap.entrySet()) {

                for (Entry<String, String> entry2 : entry1.getValue().entrySet()) {

                    String variantId = entry2.getKey();
                    String rsid = variantToRsidMap.get(variantId);

                    if (!inspectedSnp.contains(variantId) && !inspectedSnp.contains(rsid)) {

                        result.add(entry2.getValue());
                        inspectedSnp.add(variantId);

                    }

                    if (ldMatrixReader != null) {

                        ArrayList<R2> r2s = ldMatrixReader.getR2(variantId);

                        if (r2s != null && r2s.size() > 0) {

                            for (R2 r2 : r2s) {

                                if (r2.r2Value >= minR2) {

                                    String variantBId = ldMatrixReader.variantIds[r2.variantB];
                                    String variantBRsid = ldMatrixReader.getRsId(variantBId);

                                    inspectedSnp.add(variantBId);
                                    inspectedSnp.add(variantBRsid);

                                }
                            }
                        }
                    }

                    if (ensemblPopulation != null) {

                        if (rsid != null) {

                            ArrayList<ProxyCoordinates> ensemblProxies = EnsemblAPI.getProxies(
                                    rsid,
                                    ensemblPopulation,
                                    minR2,
                                    buildNumber
                            );

                            for (ProxyCoordinates proxyCoordinates : ensemblProxies) {

                                inspectedSnp.add(proxyCoordinates.proxySnp);

                            }
                        }
                    }

                    if (ldLinkPopulation != null) {

                        if (rsid != null) {

                            ArrayList<ProxyCoordinates> ldlinkProxies = LDproxy.getProxy(
                                    rsid,
                                    ldLinkPopulation,
                                    "r2",
                                    "500000",
                                    ldLinkToken
                            );

                            for (ProxyCoordinates proxyCoordinates : ldlinkProxies) {

                                if (proxyCoordinates.r2 >= minR2) {

                                    inspectedSnp.add(proxyCoordinates.proxySnp);

                                }
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
