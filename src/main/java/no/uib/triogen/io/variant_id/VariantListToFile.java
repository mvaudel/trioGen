package no.uib.triogen.io.variant_id;

import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.annotation.ensembl.EnsemblAPI;
import no.uib.triogen.model.annotation.ProxyCoordinates;
import no.uib.triogen.model.annotation.ensembl.VariantCoordinates;
import no.uib.triogen.model.annotation.ld_link.LDproxy;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * This class converts a list of variants to target file and looks for proxies
 * for missing variants.
 *
 * @author Marc Vaudel
 */
public class VariantListToFile {

    /**
     * Chomosome to rsid to position map.
     */
    private HashMap<String, HashMap<String, Integer>> rsIdCoordinatesMap = new HashMap<>(24);

    /**
     * Writes the variant file from the given variant list.
     *
     * @param source The value to put in the source column for future reference.
     * @param bgenFilePath The path to the bgen files.
     * @param variantListFile The file containing the list of variants to
     * annotate.
     * @param targetFile The file to write to.
     * @param missingFile The file where to write the SNPs that were not found.
     * @param buildNumber The build number to use.
     * @param ensemblPopulation The Ensembl population to query.
     * @param ldLinkPopulation The LDlink population to query.
     * @param ldLinkToken The token to use for LDlink.
     * @param minR2 The minimum R2 to use for proxies.
     * @param logger Logger to display feedback.
     *
     * @throws IOException Exception thrown if an error occurred while reading
     * or writing a file.
     */
    public void writeVariantFile(
            String source,
            String bgenFilePath,
            File variantListFile,
            File targetFile,
            File missingFile,
            int buildNumber,
            String ldLinkPopulation,
            String ldLinkToken,
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

                                        rsIdMap.put(variantInformation.rsid, variantInformation.position);

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

                                    ArrayList<ProxyCoordinates> ldLinkProxies = ldLinkToken != null ? new ArrayList<>(0) : LDproxy.getProxy(rsId, ldLinkPopulation, "r2", "500000", ldLinkToken);

                                    ArrayList<ProxyCoordinates> ensemblProxies = EnsemblAPI.getProxies(
                                            rsId,
                                            ensemblPopulation,
                                            minR2,
                                            buildNumber
                                    );

                                    TreeMap<Double, TreeMap<Integer, HashMap<String, ProxyCoordinates>>> proxyMap = getProxyMap(variantCoordinates, ldLinkProxies, ensemblProxies, rsIdMap);

                                    if (!proxyMap.isEmpty()) {

                                        HashMap<String, ProxyCoordinates> bestProxies = proxyMap.lastEntry().getValue().firstEntry().getValue();

                                        ProxyCoordinates bestProxy = null;

                                        for (ProxyCoordinates proxyCoordinates : bestProxies.values()) {

                                            bestProxy = proxyCoordinates;

                                            if (bestProxy.alleleMapping != null) {

                                                break;

                                            }
                                        }

                                        String description = bestProxy.alleleMapping != null
                                                ? String.join("",
                                                        "Proxy for ", rsId, " (", ensemblPopulation, " r2=", Double.toString(bestProxy.r2), " alleles=", bestProxy.alleleMapping, ")"
                                                )
                                                : String.join("",
                                                        "Proxy for ", rsId, " (", ensemblPopulation, " r2=", Double.toString(bestProxy.r2), ")"
                                                );

                                        targetWriter.writeLine(
                                                bestProxy.proxySnp,
                                                bestProxy.contig,
                                                Integer.toString(bestProxy.start),
                                                source,
                                                description
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

    /**
     * Returns a map of the proxies found in the bgen file and passing the
     * distance and r2 thresholds.
     *
     * @param variantCoordinates The coordinates of the target variant.
     * @param ldLinkProxies The proxies from LDlnk.
     * @param ensemblProxies The proxies from Ensembl.
     * @param rsIdMap A map indicating whether a given rsid is present in the
     * bgen file.
     *
     * @return A map of the proxies found in the bgen file and passing the
     * distance and r2 thresholds.
     */
    private TreeMap<Double, TreeMap<Integer, HashMap<String, ProxyCoordinates>>> getProxyMap(
            VariantCoordinates variantCoordinates,
            ArrayList<ProxyCoordinates> ldLinkProxies,
            ArrayList<ProxyCoordinates> ensemblProxies,
            HashMap<String, Integer> rsIdMap
    ) {

        TreeMap<Double, TreeMap<Integer, HashMap<String, ProxyCoordinates>>> proxyMap = new TreeMap<>();

        ArrayList<ProxyCoordinates> proxies = Stream.concat(ldLinkProxies.stream(), ensemblProxies.stream())
                .collect(Collectors.toCollection(ArrayList::new));

        for (ProxyCoordinates proxyCoordinates : proxies) {

            Integer bgenPosition = rsIdMap.get(proxyCoordinates.proxySnp);

            if (bgenPosition != null) {

                TreeMap<Integer, HashMap<String, ProxyCoordinates>> r2Map = proxyMap.get(proxyCoordinates.r2);

                if (r2Map == null) {

                    r2Map = new TreeMap<>();
                    proxyMap.put(proxyCoordinates.r2, r2Map);

                }

                int distance = Math.abs(variantCoordinates.position - bgenPosition);

                HashMap<String, ProxyCoordinates> proxiesAtDistance = r2Map.get(distance);

                if (proxiesAtDistance == null) {

                    proxiesAtDistance = new HashMap<>(1);
                    r2Map.put(distance, proxiesAtDistance);

                }

                proxiesAtDistance.put(proxyCoordinates.proxySnp, proxyCoordinates);

            }
        }

        return proxyMap;
    }
}
