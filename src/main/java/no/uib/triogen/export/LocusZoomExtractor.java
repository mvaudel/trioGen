package no.uib.triogen.export;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzReader;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import no.uib.triogen.io.ld.LdMatrixReader;
import no.uib.triogen.model.annotation.EnsemblAPI;
import no.uib.triogen.model.annotation.GeneCoordinates;

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
     * @param targetVariant The id of the variant of interest.
     * @param maxDistance The maximum distance from the variant in bp.
     * @param buildNumber The number of the build, e.g. 38 for GRCh38.
     * @param resultFile The file containing the results of the linear association.
     * @param ldFile The file containing the ld matrix.
     * @param destinationFile The file where to write the locus zoom data.
     * @param genesFile The file where to write the gene mapping.
     * 
     * @throws IOException Exception thrown if an error occurred while reading or writing a file.
     */
    public static void writeData(
            String targetPhenotype,
            String targetVariant,
            int maxDistance,
            int buildNumber,
            File resultFile,
            File ldFile,
            File destinationFile,
            File genesFile
    ) throws IOException {

        File indexFile = IoUtils.getIndexFile(resultFile);

        if (!indexFile.exists()) {

            throw new FileNotFoundException("Index file " + indexFile + " not found for " + resultFile + ".");

        }

        LdMatrixReader ldMatrixReader = new LdMatrixReader(ldFile);

        ArrayList<Integer> pValuesIndexes = new ArrayList<>();
        ArrayList<String> variables = new ArrayList<>();
        ArrayList<String> models = new ArrayList<>();

        int typedColumn = -1;
        int nColumn = -1;

        ArrayList<String> contigs = new ArrayList<>();
        ArrayList<Integer> bps = new ArrayList<>();
        ArrayList<String> variantIds = new ArrayList<>();
        ArrayList<Long> positions = new ArrayList<>();
        ArrayList<Integer> compressedLengths = new ArrayList<>();
        ArrayList<Integer> uncompressedLengths = new ArrayList<>();
        HashMap<String, Double> ldMap = ldMatrixReader.getR2(targetVariant);

        String targetContig = null;
        int targetBp = -1;

        try ( SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

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
                            writer.writeLine(resultLine.trim());

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
                        } else if (phenoName.equals(targetPhenotype)) {

                            String contig = lineSplit[0];
                            int bp = Integer.parseInt(lineSplit[1]);

                            if (variantId.equals(targetVariant)) {

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

                                targetContig = contig;
                                targetBp = bp;

                                for (int i = 0; i < contigs.size(); i++) {

                                    contig = contigs.get(i);

                                    if (contig.equals(targetContig)) {

                                        bp = bps.get(i);

                                        if (Math.abs(bp - targetBp) <= maxDistance) {

                                            variantId = variantIds.get(i);
                                            double ld = ldMap.getOrDefault(variantId, 0.0);

                                            long previousPosition = positions.get(i);
                                            int previousCompressedLength = compressedLengths.get(i);
                                            int previousUncompressedLength = uncompressedLengths.get(i);

                                            String resultLine = gzReader.read(previousPosition, previousCompressedLength, previousUncompressedLength);

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
                                            targetContig,
                                            Integer.toString(targetBp),
                                            targetVariant,
                                            typed,
                                            n,
                                            Double.toString(1.0),
                                            model,
                                            variable,
                                            Double.toString(pValue)
                                    );
                                }

                            } else if (targetContig == null) {

                                contigs.add(contig);
                                bps.add(bp);
                                variantIds.add(variantId);
                                positions.add(position);
                                compressedLengths.add(compressedLength);
                                uncompressedLengths.add(uncompressedLength);

                            } else {

                                if (Math.abs(bp - targetBp) <= maxDistance) {

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
        }

        if (genesFile != null && targetContig != null) {

            try ( SimpleFileWriter writer = new SimpleFileWriter(genesFile, true)) {

                String ensemblVersion = EnsemblAPI.getEnsemblVersion(buildNumber);

                writer.writeLine("# Ensembl version: " + ensemblVersion);
                writer.writeLine("biotype", "name", "start", "end");

                if (maxDistance > 2.5e6) {

                    throw new IllegalArgumentException("Maximal window size for gene coordinates is 5e6.");

                }

                ArrayList<GeneCoordinates> geneCoordinatesList = EnsemblAPI.getGeneCoordinates(
                        targetContig, 
                        targetBp - maxDistance, 
                        targetBp + maxDistance,
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
