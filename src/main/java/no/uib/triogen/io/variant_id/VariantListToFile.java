package no.uib.triogen.io.variant_id;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.annotation.EnsemblAPI;
import no.uib.triogen.model.annotation.ProxyCoordinates;
import no.uib.triogen.model.annotation.VariantCoordinates;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * This class converts a list of variants to target file and looks for proxies
 * for missing variants.
 *
 * @author Marc Vaudel
 */
public class VariantListToFile {

    private HashMap<String, HashMap<String, Integer>> rsIdCoordinatesMap = new HashMap<>(24);

    public void writeVariantFile(
            String source,
            String bgenFilePath,
            File variantListFile,
            File targetFile,
            File missingFile,
            int buildNumber,
            String ensemblPopulation,
            double minR2,
            SimpleCliLogger logger
    ) throws IOException {

        logger.logMessage("Mapping variants from " + variantListFile + ".");

        long start = Instant.now().getEpochSecond();

        int nVariants = 0;
        int nFound = 0;

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(variantListFile)) {

            try (SimpleFileWriter targetWriter = new SimpleFileWriter(targetFile, false)) {

                targetWriter.writeLine("snp", "chr", "position", "source", "comment");

                try (SimpleFileWriter missingWriter = new SimpleFileWriter(missingFile, false)) {

                    String rsId;
                    while ((rsId = reader.readLine()) != null) {

                        logger.logMessage("Mapping " + rsId + ".");

                        nVariants++;

                        if (!rsId.startsWith("rs")) {

                            throw new IllegalArgumentException("Variant id " + rsId + " not supported, please provide an rs id.");

                        }

                        ArrayList<VariantCoordinates> variantCoordinatesArray = EnsemblAPI.getVariantCoordinates(rsId, buildNumber);

                        if (variantCoordinatesArray.size() > 1) {

                            logger.logMessage("More than one genomic coordinate (" + variantCoordinatesArray.size() + ") found for variant " + rsId + ", all will be included.");

                        } else if (variantCoordinatesArray.isEmpty()) {

                            missingWriter.writeLine(rsId);

                        } else {

                            for (VariantCoordinates variantCoordinates : variantCoordinatesArray) {

                                String chr = variantCoordinates.contig;

                                File bgenFile = new File(bgenFilePath.replace("{chr}", chr));

                                if (!bgenFile.exists() && variantCoordinates.contig.equals("X")) {

                                    chr = "23";
                                    bgenFile = new File(bgenFilePath.replace("{chr}", chr));

                                } else if (!bgenFile.exists() && variantCoordinates.contig.equals("Y")) {

                                    chr = "24";
                                    bgenFile = new File(bgenFilePath.replace("{chr}", chr));

                                }

                                if (!bgenFile.exists()) {

                                    throw new IllegalArgumentException("Bgen file not found for chromosome " + chr + ".");

                                }
                                logger.logMessage("Parsing bgen file for chromosome " + chr + ".");

                                HashMap<String, Integer> rsIdMap = rsIdCoordinatesMap.get(chr);

                                if (rsIdMap == null) {

                                    BgenIndex bgenIndex = BgenIndex.getBgenIndex(bgenFile);
                                    rsIdMap = new HashMap<>(bgenIndex.variantIdArray.length);
                                    rsIdCoordinatesMap.put(chr, rsIdMap);

                                    for (VariantInformation variantInformation : bgenIndex.variantInformationArray) {

                                        rsIdMap.put(variantInformation.rsId, variantInformation.position);

                                    }
                                }

                                Integer position = rsIdMap.get(rsId);

                                if (position != null) {

                                    targetWriter.writeLine(
                                            rsId,
                                            chr,
                                            position.toString(),
                                            source,
                                            ""
                                    );

                                    logger.logMessage("Found " + rsId + " in gben files (chr " + chr + ", pos " + position + ").");

                                    nFound++;

                                } else {

                                    ArrayList<ProxyCoordinates> proxies = EnsemblAPI.getProxies(
                                            rsId,
                                            ensemblPopulation,
                                            minR2,
                                            buildNumber
                                    );

                                    TreeMap<Double, TreeMap<Integer, ArrayList<ProxyCoordinates>>> proxyMap = new TreeMap<>();

                                    for (ProxyCoordinates proxyCoordinates : proxies) {

                                        Integer bgenPosition = rsIdMap.get(proxyCoordinates.proxySnp);

                                        if (bgenPosition != null) {

                                            TreeMap<Integer, ArrayList<ProxyCoordinates>> r2Map = proxyMap.get(proxyCoordinates.r2);

                                            if (r2Map == null) {

                                                r2Map = new TreeMap<>();
                                                proxyMap.put(proxyCoordinates.r2, r2Map);

                                            }

                                            int distance = Math.abs(variantCoordinates.position - position);

                                            ArrayList<ProxyCoordinates> proxiesAtDistance = r2Map.get(distance);

                                            if (proxiesAtDistance == null) {

                                                proxiesAtDistance = new ArrayList<>(1);
                                                r2Map.put(distance, proxiesAtDistance);

                                            }

                                            proxiesAtDistance.add(proxyCoordinates);

                                        }
                                    }

                                    if (!proxyMap.isEmpty()) {

                                        ProxyCoordinates bestProxy = proxyMap.lastEntry().getValue().firstEntry().getValue().get(0);

                                        targetWriter.writeLine(
                                                bestProxy.proxySnp,
                                                bestProxy.contig,
                                                Integer.toString(bestProxy.start),
                                                source,
                                                String.join("",
                                                        "Proxy for ", rsId, " (", ensemblPopulation, " r2=", Double.toString(bestProxy.r2), ")"
                                                )
                                        );

                                        logger.logMessage("Proxy found for " + rsId + " (r2=" + bestProxy.r2 + ").");

                                        nFound++;

                                    } else {

                                        missingWriter.writeLine(rsId);

                                        logger.logMessage("No proxy found for " + rsId + " (chr " + chr + ", pos " + position + ").");

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        long end = Instant.now().getEpochSecond();
        long duration = end - start;

        logger.logMessage("Mapping variants done (" + nFound + " variants mapped from " + nVariants + " in " + duration + " seconds)");

    }
}
