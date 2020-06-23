package no.uib.triogen.export;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzReader;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.annotation.EnsemblAPI;
import no.uib.triogen.model.annotation.GeneCoordinates;
import no.uib.triogen.model.trio_genotypes.VariantList;

/**
 * This class extracts the data necessary to build locus zoom plots.
 *
 * @author Marc Vaudel
 */
public class LocusZoomExtractor {

    /**
     * Writes the data needed for the given locus zoom.
     *
     * @param targetPhenotype The phenotype of interest.
     * @param variantList The ids of the variant of interest.
     * @param maxDistance The maximum distance from the variant in bp.
     * @param buildNumber The number of the build, e.g. 38 for GRCh38.
     * @param resultFile The file containing the results of the linear
     * association.
     * @param ldFile The file containing the ld matrix.
     * @param destinationFileStem The file where to write the locus zoom data.
     * @param genesFileStem The file where to write the gene mapping.
     * @param logger The logger to use.
     *
     * @throws IOException Exception thrown if an error occurred while reading
     * or writing a file.
     */
    public static void writeData(
            String targetPhenotype,
            VariantList variantList,
            int maxDistance,
            int buildNumber,
            File resultFile,
            File ldFile,
            String destinationFileStem,
            String genesFileStem,
            SimpleCliLogger logger
    ) throws IOException {

        logger.logMessage("Getting LD values between the target variants and nearby variants.");

        HashMap<String, HashMap<String, Double>> ldMaps = getLdMaps(variantList, ldFile);

        logger.logMessage("Iterating association results");

        File indexFile = IoUtils.getIndexFile(resultFile);

        if (!indexFile.exists()) {

            throw new FileNotFoundException("Index file " + indexFile + " not found for " + resultFile + ".");

        }

        ArrayList<Integer> pValuesIndexes = new ArrayList<>();
        ArrayList<String> variables = new ArrayList<>();
        ArrayList<String> models = new ArrayList<>();

        int typedColumn = -1;
        int nColumn = -1;

        HashSet<Integer> variantsFound = new HashSet<>();

        HashMap<String, HashMap<String, SimpleFileWriter>> fileWriters = new HashMap<>(variantList.variantId.length);

        ArrayList<String> comments = new ArrayList<>();

        long position = IndexedGzWriter.HEADER_LENGTH;

        try ( IndexedGzReader gzReader = new IndexedGzReader(resultFile)) {

            try ( SimpleFileReader indexReader = SimpleFileReader.getFileReader(indexFile)) {

                String indexLine = indexReader.readLine();

                while ((indexLine = indexReader.readLine()) != null) {

                    String[] lineSplit = indexLine.split(IoUtils.SEPARATOR);
                    String variantId = lineSplit[2];
                    String phenoName = lineSplit[3];
                    int compressedLength = Integer.parseInt(lineSplit[4]);
                    int uncompressedLength = Integer.parseInt(lineSplit[5]);

                    if (phenoName.equals("Comment")) {

                        String resultLine = gzReader.read(position, compressedLength, uncompressedLength);
                        comments.add(resultLine);

                    } else if (phenoName.equals("Header")) {

                        String resultLine = gzReader.read(position, compressedLength, uncompressedLength);
                        lineSplit = resultLine
                                .trim()
                                .split(IoUtils.SEPARATOR);

                        for (int i = 0; i < lineSplit.length; i++) {

                            if (lineSplit[i].endsWith(".p")) {

                                String temp = lineSplit[i];
                                temp = temp.substring(0, temp.length() - 2);
                                String[] subSplit = temp.split("\\.");

                                String model = subSplit[0];
                                String variable = subSplit[1];

                                pValuesIndexes.add(i);
                                models.add(model);
                                variables.add(variable);

                            } else if (lineSplit[i].equals("typed")) {

                                typedColumn = i;

                            } else if (lineSplit[i].equals("n")) {

                                nColumn = i;

                            }
                        }
                    } else if (targetPhenotype == null || phenoName.equals(targetPhenotype)) {

                        String contig = lineSplit[0];
                        int bp = Integer.parseInt(lineSplit[1]);

                        for (int variantI = 0; variantI < variantList.variantId.length; variantI++) {

                            if (contig.equals(variantList.chromosome[variantI]) && bp >= variantList.start[variantI] - maxDistance && bp <= variantList.end[variantI] + maxDistance) {

                                String targetVariantId = variantList.variantId[variantI];

                                HashMap<String, SimpleFileWriter> variantWritersMap = fileWriters.get(targetVariantId);

                                if (variantWritersMap == null) {

                                    variantWritersMap = new HashMap<>(1);
                                    fileWriters.put(targetVariantId, variantWritersMap);

                                }

                                SimpleFileWriter writer = variantWritersMap.get(phenoName);

                                if (writer == null) {

                                    variantsFound.add(variantI);

                                    File destinationFile = new File(destinationFileStem + "_" + variantId + "_" + phenoName + "_LocusZoomData.gz");
                                    writer = new SimpleFileWriter(destinationFile, true);

                                    for (String comment : comments) {

                                        writer.writeLine(comment);

                                    }

                                    writer.writeLine(
                                            "contig",
                                            "position",
                                            "variantId",
                                            "typed",
                                            "n",
                                            "ld",
                                            "model",
                                            "variable",
                                            "p"
                                    );

                                    variantWritersMap.put(phenoName, writer);

                                }

                                HashMap<String, Double> ldMap = ldMaps.get(targetVariantId);
                                double ld = ldMap.getOrDefault(variantId, 0.0);

                                String resultLine = gzReader.read(position, compressedLength, uncompressedLength);

                                lineSplit = resultLine
                                        .trim()
                                        .split(IoUtils.SEPARATOR);

                                String typed = lineSplit[typedColumn];
                                String n = lineSplit[nColumn];

                                for (int j = 0; j < pValuesIndexes.size(); j++) {

                                    double pValue = Double.parseDouble(lineSplit[pValuesIndexes.get(j)]);
                                    String model = models.get(j);
                                    String variable = variables.get(j);

                                    writer.writeLine(
                                            contig,
                                            Integer.toString(bp),
                                            variantId,
                                            typed,
                                            n,
                                            Double.toString(ld),
                                            model,
                                            variable,
                                            Double.toString(pValue)
                                    );
                                }
                            }
                        }
                    }

                    position += compressedLength;

                }
            }
        }

        fileWriters.values().stream()
                .flatMap(
                        variantMap -> variantMap.values().stream()
                )
                .forEach(
                        writer -> writer.close()
                );

        if (genesFileStem != null && !variantsFound.isEmpty()) {

            logger.logMessage("Getting gene mapping");

            for (int variantI : variantsFound) {

                String targetContig = variantList.chromosome[variantI];
                int targetBpStart = variantList.start[variantI];
                int targetBpEnd = variantList.end[variantI];
                String targetVariantId = variantList.variantId[variantI];

                File geneFile = new File(genesFileStem + "_" + targetVariantId + "_LocusZoomGenes.gz");

                try ( SimpleFileWriter writer = new SimpleFileWriter(geneFile, true)) {

                    String ensemblVersion = EnsemblAPI.getEnsemblVersion(buildNumber);

                    writer.writeLine("# Ensembl version: " + ensemblVersion);
                    writer.writeLine("biotype", "name", "start", "end");

                    if (maxDistance > 2.5e6) {

                        throw new IllegalArgumentException("Maximal window size for gene coordinates is 5e6.");

                    }

                    ArrayList<GeneCoordinates> geneCoordinatesList = EnsemblAPI.getGeneCoordinates(
                            targetContig,
                            targetBpStart - maxDistance,
                            targetBpEnd + maxDistance,
                            buildNumber
                    );

                    for (GeneCoordinates geneCoordinates : geneCoordinatesList) {

                        writer.writeLine(
                                geneCoordinates.biotype,
                                geneCoordinates.name,
                                Integer.toString(geneCoordinates.start),
                                Integer.toString(geneCoordinates.end)
                        );
                    }
                }
            }
        }
    }

    /**
     * Returns the variants in LD with all variants in the variant list and the associated r2.
     * 
     * @param variantList The ids of the variant of interest.
     * @param ldFile The file containing the ld matrix.
     * 
     * @return The variants in LD with all variants in the variant list and the associated r2.
     * 
     * @throws IOException Exception thrown if an error occurred while retrieving the ld.
     */
    private static HashMap<String, HashMap<String, Double>> getLdMaps(
            VariantList variantList,
            File ldFile
    ) throws IOException {

        LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);

        HashMap<String, HashMap<String, Double>> result = new HashMap<>(variantList.variantId.length);

        for (int i = 0; i < variantList.variantId.length; i++) {

            String variantId = variantList.variantId[i];

            HashMap<String, Double> ldMap = ldMatrixReader.getR2(variantId);

            if (ldMap == null) {

                ldMap = new HashMap<>(0);

            }

            result.put(variantId, ldMap);

        }

        return result;

    }
}
