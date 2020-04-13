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

/**
 * This class extracts the data necessary to build locus zoom plots.
 *
 * @author Marc Vaudel
 */
public class LocusZoomExtractor {

    public static void writeData(
            String phenotype,
            String targetVariant,
            int maxDistance,
            File resultFile,
            File ldFile,
            File destinationFile
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
        HashMap<String, Double> ldMap = null;

        String targetContig = null;
        int targetBp = -1;

        try ( SimpleFileWriter writer = new SimpleFileWriter(destinationFile, true)) {

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
                            writer.writeLine(resultLine);

                        } else if (phenoName.equals("Header")) {

                            String resultLine = gzReader.read(position, compressedLength, uncompressedLength);
                            lineSplit = resultLine
                                    .trim()
                                    .split(IoUtils.SEPARATOR);

                            for (int i = 0; i < lineSplit.length; i++) {

                                if (lineSplit[i].endsWith("_p")) {

                                    String temp = lineSplit[i];
                                    temp = temp.substring(0, temp.length() - 2);
                                    int sep = temp.lastIndexOf('_');

                                    String model = temp.substring(0, sep);
                                    String variable = temp.substring(sep + 1);

                                    pValuesIndexes.add(i);
                                    models.add(model);
                                    variables.add(variable);

                                } else if (lineSplit[i].equals("typed")) {

                                    typedColumn = i;

                                } else if (lineSplit[i].equals("n")) {

                                    nColumn = i;

                                }
                            }
                        } else {

                            String contig = lineSplit[0];
                            int bp = Integer.parseInt(lineSplit[1]);

                            if (variantId.equals(targetVariant)) {

                                targetContig = contig;
                                targetBp = bp;

                                ldMap = ldMatrixReader.getR2(targetVariant);

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
                    }
                }
            }
        }
    }
}
