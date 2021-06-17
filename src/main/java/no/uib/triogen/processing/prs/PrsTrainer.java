package no.uib.triogen.processing.prs;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.SEPARATOR;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.ld.R2;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.utils.Utils;

/**
 * This class trains a PRS on trios.
 *
 * @author Marc Vaudel
 */
public class PrsTrainer {

    /**
     * Wildcard for the variable name in the summary stats columns.
     */
    public static final String VARIABLE_WILDCARD = "{variable}";
    /**
     * The ld matrix file path.
     */
    private final String ldMatrixFilePath;
    /**
     * Map of the LD matrix readers.
     */
    private final HashMap<String, LdMatrixReader> ldMatrixReadersMap = new HashMap<>(24);
    /**
     * The file to export the score details to.
     */
    private final File destinationFile;
    /**
     * The file containing the training summary stats.
     */
    private final File trainingFile;
    /**
     * The column containing the snp id.
     */
    private final String snpIdColumn;
    /**
     * The column containing the chromosome.
     */
    private final String chrColumn;
    /**
     * The column containing the position.
     */
    private final String posColumn;
    /**
     * The column containing the reference allele.
     */
    private final String refColumn;
    /**
     * The column containing the effect allele.
     */
    private final String eaColumn;
    /**
     * The pattern to use to find the effect size column for each variable in
     * the model.
     */
    private final String betaColumnPattern;
    /**
     * The pattern to use to find the standard error column for each variable in
     * the model.
     */
    private final String seColumnPattern;
    /**
     * The pattern to use to find the p-value column for each variable in the
     * model.
     */
    private final String pValueColumnPattern;
    /**
     * The trio model to use.
     */
    private final Model model;
    /**
     * The ordered names of the variables used in the model.
     */
    private final String[] variableNames;
    /**
     * The LD threshold used for pruning a locus.
     */
    private final double ldLocusThreshold;
    /**
     * The LD threshold used around a top hit.
     */
    private final double ldTopHitThreshold;
    /**
     * The minimal number of SNPs to consider in a locus.
     */
    private final int nSnpPerLocusThreshold;
    /**
     * The logger.
     */
    private final SimpleCliLogger logger;

    /**
     * Constructor.
     *
     * @param trainingFile The file containing the training data.
     * @param ldMatrixFilePath The LD matrix file path.
     * @param destinationFile The file to export the result to.
     * @param snpIdColumn The column containing SNP ids.
     * @param chrColumn The column containing the chromosome.
     * @param posColumn The column containing the position.
     * @param refColumn The column containing the reference allele.
     * @param eaColumn The name of the effect allele column.
     * @param model The trio model to use.
     * @param variableNames The names of the variables to include.
     * @param betaColumnPattern The pattern to use to find the effect size
     * column for each variable in the model.
     * @param seColumnPattern The pattern to use to find the standard error
     * column for each variable in the model.
     * @param pValueColumnPattern The pattern to use to find the effect size
     * column for each variable in the model.
     * @param ldLocusThreshold The LD threshold used for pruning.
     * @param ldTopHitThreshold The LD threshold used for top hits.
     * @param nSnpPerLocusThreshold The number of variants to consider per
     * locus.
     * @param logger The logger.
     */
    public PrsTrainer(
            File trainingFile,
            String ldMatrixFilePath,
            File destinationFile,
            String snpIdColumn,
            String chrColumn,
            String posColumn,
            String refColumn,
            String eaColumn,
            String betaColumnPattern,
            String seColumnPattern,
            String pValueColumnPattern,
            Model model,
            String[] variableNames,
            int nSnpPerLocusThreshold,
            double ldLocusThreshold,
            double ldTopHitThreshold,
            SimpleCliLogger logger
    ) {

        this.trainingFile = trainingFile;
        this.ldMatrixFilePath = ldMatrixFilePath;
        this.destinationFile = destinationFile;
        this.snpIdColumn = snpIdColumn;
        this.chrColumn = chrColumn;
        this.posColumn = posColumn;
        this.refColumn = refColumn;
        this.eaColumn = eaColumn;
        this.betaColumnPattern = betaColumnPattern;
        this.seColumnPattern = seColumnPattern;
        this.pValueColumnPattern = pValueColumnPattern;
        this.model = model;
        this.variableNames = variableNames;
        this.nSnpPerLocusThreshold = nSnpPerLocusThreshold;
        this.ldLocusThreshold = ldLocusThreshold;
        this.ldTopHitThreshold = ldTopHitThreshold;
        this.logger = logger;

    }

    public void run() throws IOException {

        logger.logMessage("Parsing training data from " + trainingFile.getAbsolutePath());

        long start = Instant.now().getEpochSecond();

        TrainingData trainingData = parseTrainingData();

        TreeMap<Double, TreeSet<String>> pValueToSnpMap = new TreeMap<>();

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Parsing training data from " + trainingFile + " done (" + duration + " seconds)");

        logger.logMessage("Weighted pruning");

        start = Instant.now().getEpochSecond();

        pruneAndExport(trainingData);

        end = Instant.now().getEpochSecond();
        duration = end - start;

        logger.logMessage("Parsing training data from " + trainingFile + " done (" + duration + " seconds)");

    }

    private void pruneAndExport(TrainingData trainingData) {

        try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

            StringBuilder line = new StringBuilder();

            line.append("variantid")
                    .append(IoUtils.SEPARATOR)
                    .append("rsid")
                    .append(IoUtils.SEPARATOR)
                    .append("chr")
                    .append(IoUtils.SEPARATOR)
                    .append("pos")
                    .append(IoUtils.SEPARATOR)
                    .append("ref_allele")
                    .append(IoUtils.SEPARATOR)
                    .append("effect_allele")
                    .append(IoUtils.SEPARATOR)
                    .append("weight");

            boolean singlePValue = !pValueColumnPattern.contains(VARIABLE_WILDCARD);

            for (int j = 0; j < variableNames.length; j++) {

                String variable = variableNames[j];

                String betaColumn = betaColumnPattern
                        .replace(VARIABLE_WILDCARD, variable);

                String seValueColumn = seColumnPattern
                        .replace(VARIABLE_WILDCARD, variable);

                line.append(IoUtils.SEPARATOR)
                        .append(betaColumn)
                        .append(IoUtils.SEPARATOR)
                        .append(seValueColumn);

                if (!singlePValue) {

                    String pValueColumn = pValueColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    line.append(IoUtils.SEPARATOR)
                            .append(pValueColumn);

                }
            }

            if (singlePValue) {

                line.append(IoUtils.SEPARATOR)
                        .append(pValueColumnPattern);

            }

            writer.writeLine(line.toString());

            HashSet<String> processedVariants = new HashSet<>();

            for (Entry<Double, TreeSet<String>> pEntry : trainingData.pValueToVariantMap.entrySet()) {

                for (String variantId : pEntry.getValue()) {

                    if (!processedVariants.contains(variantId)) {

                        String[] variantDetails = trainingData.variantToDetailsMap.get(variantId);

                        LdMatrixReader ldMatrixReader = getLdMatrixReader(variantDetails[0]);

                        ArrayList<R2> r2s = ldMatrixReader.getR2(variantId);

                        ArrayList<R2> r2InLocus = new ArrayList<>();
                        HashSet<String> idsInLocus = new HashSet<>();

                        for (R2 r2 : r2s) {

                            if (r2.r2Value >= ldLocusThreshold) {

                                r2InLocus.add(r2);
                                idsInLocus.add(r2.getVariantBId());
                                idsInLocus.add(r2.getVariantBRsid());

                            }
                        }

                        if (r2InLocus.size() >= nSnpPerLocusThreshold) {

                            String[] bestSnps = new String[variantDetails.length];
                            String[] bestRsids = new String[variantDetails.length];
                            double[] bestPs = new double[variantDetails.length];

                            for (int variableI = 0; variableI < variantDetails.length; variableI++) {

                                bestPs[variableI] = 1.1;

                                HashMap<String, double[]> variableResult = trainingData.variantToSummaryStats[variableI];

                                for (R2 r2 : r2InLocus) {

                                    double[] summaryStats = variableResult.get(r2.getVariantBId());

                                    if (summaryStats != null && summaryStats[2] < bestPs[variableI]) {

                                        bestSnps[variableI] = r2.getVariantBId();
                                        bestRsids[variableI] = r2.getVariantBRsid();
                                        bestPs[variableI] = summaryStats[2];

                                    }

                                    summaryStats = variableResult.get(r2.getVariantBRsid());

                                    if (summaryStats != null && summaryStats[2] < bestPs[variableI]) {

                                        bestSnps[variableI] = r2.getVariantBRsid();
                                        bestPs[variableI] = summaryStats[2];

                                    }
                                }
                            }

                            HashMap<String, String> topHits = new HashMap<>();
                            HashMap<String, String> hitsIds = new HashMap<>();

                            for (int variableI = 0; variableI < variantDetails.length; variableI++) {

                                String topHitId = bestSnps[variableI];
                                String topHitRsId = bestRsids[variableI];

                                topHits.put(topHitId, topHitRsId);

                                r2s = ldMatrixReader.getR2(topHitId);

                                for (R2 r2 : r2s) {

                                    if (r2.r2Value >= ldTopHitThreshold) {

                                        if (trainingData.variantToDetailsMap.containsKey(r2.getVariantBId())) {

                                            topHits.put(r2.getVariantBId(), r2.getVariantBRsid());
                                            hitsIds.put(r2.getVariantBId(), r2.getVariantBId());

                                        } else if (trainingData.variantToDetailsMap.containsKey(r2.getVariantBRsid())) {

                                            topHits.put(r2.getVariantBId(), r2.getVariantBRsid());
                                            hitsIds.put(r2.getVariantBId(), r2.getVariantBId());

                                        }
                                    }
                                }
                            }

                            for (Entry<String, String> entry : topHits.entrySet()) {

                                String hitId = entry.getKey();
                                String rsId = entry.getValue();

                                String id = hitsIds.get(hitId);
                                String[] hitDetails = trainingData.variantToDetailsMap.get(id);

                                double weight = 1.0 / topHits.size();

                                line = new StringBuilder();

                                line.append(hitId)
                                        .append(IoUtils.SEPARATOR)
                                        .append(rsId)
                                        .append(IoUtils.SEPARATOR)
                                        .append(hitDetails[0])
                                        .append(IoUtils.SEPARATOR)
                                        .append(hitDetails[1])
                                        .append(IoUtils.SEPARATOR)
                                        .append(hitDetails[2])
                                        .append(IoUtils.SEPARATOR)
                                        .append(hitDetails[3])
                                        .append(IoUtils.SEPARATOR)
                                        .append(weight);

                                double[] summaryStats = null;

                                for (int j = 0; j < variableNames.length; j++) {

                                    summaryStats = trainingData.variantToSummaryStats[j].get(id);

                                    line.append(IoUtils.SEPARATOR)
                                            .append(summaryStats[0])
                                            .append(IoUtils.SEPARATOR)
                                            .append(summaryStats[1]);

                                    if (!singlePValue) {

                                        line.append(IoUtils.SEPARATOR)
                                                .append(summaryStats[0]);

                                    }
                                }

                                if (singlePValue) {

                                    line.append(IoUtils.SEPARATOR)
                                            .append(summaryStats[2]);

                                }

                                writer.writeLine(line.toString());

                            }
                        }

                        processedVariants.addAll(idsInLocus);

                    }
                }
            }
        }
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

        LdMatrixReader ldMatrixReader = ldMatrixReadersMap.get(contig);

        if (ldMatrixReader == null) {

            String contigLdMatrixFilePath = ldMatrixFilePath.replace(Utils.CONTIG_WILDCARD, contig);

            File contigLdFile = new File(contigLdMatrixFilePath);

            if (!contigLdFile.exists()) {

                throw new IllegalArgumentException("No LD matrix found for chromosome " + contig + " at " + contigLdMatrixFilePath + ".");

            }

            try {

                ldMatrixReader = new LdMatrixReader(contigLdFile);

            } catch (Exception e) {

                throw new RuntimeException(e);

            }

        }

        return ldMatrixReader;

    }

    /**
     * For each model variable, parses the training data in a map, variant to
     * effect and p.
     *
     * @return The training data in a map, variant to effect and p, per
     * variable.
     */
    private TrainingData parseTrainingData() {

        HashMap<String, double[]>[] variantToSummaryStats = new HashMap[variableNames.length];
        TreeMap<Double, TreeSet<String>> pValueToVariantMap = new TreeMap<>();
        HashMap<String, String[]> variantToDetailsMap = new HashMap<>();

        for (int j = 0; j < variableNames.length; j++) {

            variantToSummaryStats[j] = new HashMap<>();

        }

        int snpIdIndex = -1;
        int chrIndex = -1;
        int posIndex = -1;
        int refAlleleIndex = -1;
        int testedAlleleIndex = -1;
        int[] betaColumnIndexes = new int[variableNames.length];
        int[] seColumnIndexes = new int[variableNames.length];
        int[] pColumnIndexes = new int[variableNames.length];

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(trainingFile)) {

            String line = reader.readLine();

            String[] lineSplit = line.split(SEPARATOR);

            for (int i = 0; i < lineSplit.length; i++) {

                if (lineSplit[i].equals(snpIdColumn)) {

                    snpIdIndex = i;

                }

                if (lineSplit[i].equals(chrColumn)) {

                    chrIndex = i;

                }

                if (lineSplit[i].equals(posColumn)) {

                    posIndex = i;

                }

                if (lineSplit[i].equals(refColumn)) {

                    refAlleleIndex = i;

                }

                if (lineSplit[i].equals(eaColumn)) {

                    testedAlleleIndex = i;

                }

                for (int j = 0; j < variableNames.length; j++) {

                    String variable = variableNames[j];

                    String betaColumn = betaColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(betaColumn)) {

                        betaColumnIndexes[j] = i;

                    }

                    String seValueColumn = seColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(seValueColumn)) {

                        seColumnIndexes[j] = i;

                    }

                    String pValueColumn = pValueColumnPattern
                            .replace(VARIABLE_WILDCARD, variable);

                    if (lineSplit[i].equals(pValueColumn)) {

                        pColumnIndexes[j] = i;

                    }
                }
            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(SEPARATOR);

                String snpId = lineSplit[snpIdIndex];
                String chr = lineSplit[chrIndex];
                String pos = lineSplit[posIndex];
                String ref = lineSplit[refAlleleIndex];
                String alt = lineSplit[testedAlleleIndex];

                for (int j = 0; j < variableNames.length; j++) {

                    String betaString = lineSplit[betaColumnIndexes[j]];
                    String seString = lineSplit[seColumnIndexes[j]];
                    String pString = lineSplit[pColumnIndexes[j]];

                    double beta = Double.parseDouble(betaString);
                    double se = Double.parseDouble(seString);
                    double p = Double.parseDouble(pString);

                    variantToSummaryStats[j].put(snpId, new double[]{beta, se, p});

                    TreeSet<String> variantsAtPvalue = pValueToVariantMap.get(p);

                    if (variantsAtPvalue == null) {

                        variantsAtPvalue = new TreeSet<>();
                        pValueToVariantMap.put(p, variantsAtPvalue);

                    }

                    variantsAtPvalue.add(snpId);

                    variantToDetailsMap.put(snpId, new String[]{chr, pos, ref, alt});

                }
            }
        }

        return new TrainingData(
                variantToSummaryStats,
                pValueToVariantMap,
                variantToDetailsMap
        );

    }

    /**
     * Place holder for the results of the GWAS training data.
     */
    private class TrainingData {

        public final HashMap<String, double[]>[] variantToSummaryStats;
        public final TreeMap<Double, TreeSet<String>> pValueToVariantMap;
        public final HashMap<String, String[]> variantToDetailsMap;

        public TrainingData(
                HashMap<String, double[]>[] variantToSummaryStats,
                TreeMap<Double, TreeSet<String>> pValueToVariantMap,
                HashMap<String, String[]> variantToDetailsMap
        ) {

            this.variantToSummaryStats = variantToSummaryStats;
            this.pValueToVariantMap = pValueToVariantMap;
            this.variantToDetailsMap = variantToDetailsMap;

        }
    }
}
