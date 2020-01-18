package no.uib.triogen.cmd.results;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.lineSeparator;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzReader;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;

/**
 * Runs multiple linear models for the association with phenotypes.
 *
 * @author Marc Vaudel
 */
public class Extract {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        if (args.length == 0
                || args.length == 1 && args[0].equals("-h")
                || args.length == 1 && args[0].equals("--help")) {

            printHelp();
            return;

        }

        if (args.length == 1 && args[0].equals("-v")
                || args.length == 1 && args[0].equals("--version")) {

            System.out.println(TrioGen.getVersion());

            return;

        }

        try {

            Options lOptions = new Options();
            ExtractOptions.createOptionsCLI(lOptions);
            CommandLineParser parser = new DefaultParser();
            CommandLine commandLine = parser.parse(lOptions, args);

            ExtractOptionsBean bean = new ExtractOptionsBean(commandLine);

            run(bean);

        } catch (Throwable e) {

            e.printStackTrace();
        }
    }

    /**
     * Runs the command.
     *
     * @param bean the bean of command line parameters
     *
     * @throws IOException Exception thrown if an I/O error occurs.
     */
    private static void run(
            ExtractOptionsBean bean
    ) throws IOException {

        HashMap<String, IndexedGzWriter> gzWriters = new HashMap<>();
        HashMap<String, SimpleFileWriter> indexWriters = new HashMap<>();

        ArrayList<String> headerComments = new ArrayList<>(2);
        String headerLine = null;

        ArrayList<String> columns = new ArrayList<>();
        HashMap<String, Integer> columnIndexes = new HashMap<>();

        File resultFile = bean.inputFile;
        long position = IndexedGzWriter.HEADER_LENGTH;

        try (IndexedGzReader gzReader = new IndexedGzReader(resultFile)) {

            File indexFile = IoUtils.getIndexFile(resultFile);

            if (!indexFile.exists()) {

                throw new RuntimeException(
                        new FileNotFoundException("Index file " + indexFile + " not found for " + resultFile + ".")
                );
            }

            try (SimpleFileReader indexReader = SimpleFileReader.getFileReader(indexFile)) {

                HashSet<String> variantIds = bean.variantIds == null
                        ? new HashSet<>(0)
                        : Arrays.stream(bean.variantIds)
                                .collect(Collectors.toCollection(HashSet::new));
                HashSet<String> phenoNames = bean.phenoNames == null
                        ? new HashSet<>(0)
                        : Arrays.stream(bean.phenoNames)
                                .collect(Collectors.toCollection(HashSet::new));

                String indexLine = indexReader.readLine();

                while ((indexLine = indexReader.readLine()) != null) {

                    String[] lineSplit = indexLine.split(IoUtils.separator);

                    String variantId = lineSplit[0];
                    String phenoName = lineSplit[1];
                    int compressedLength = Integer.parseInt(lineSplit[2]);
                    int uncompressedLength = Integer.parseInt(lineSplit[3]);

                    position += compressedLength;

                    if (phenoName.equals("Comment")) {

                        String resultLine = gzReader.read(position, compressedLength, uncompressedLength);
                        headerComments.add(resultLine);

                    } else if (phenoName.equals("Header")) {

                        String resultLine = gzReader.read(position, compressedLength, uncompressedLength);

                        String[] headerSplit = resultLine
                                .split(IoUtils.separator);

                        if (bean.columns != null) {

                            HashSet<String> valuesSet = bean.columns == null ? new HashSet<>(0)
                                    : Arrays.stream(bean.columns)
                                            .collect(
                                                    Collectors.toCollection(HashSet::new)
                                            );

                            HashSet<String> headerSet = Arrays.stream(headerSplit)
                                    .collect(
                                            Collectors.toCollection(HashSet::new)
                                    );

                            String missingColumns = valuesSet.stream()
                                    .filter(
                                            column -> !headerSet.contains(column)
                                    )
                                    .collect(
                                            Collectors.joining(",")
                                    );

                            if (missingColumns.length() > 0) {

                                throw new IllegalArgumentException("Column not found: " + missingColumns + ".");

                            }

                            for (int i = 0; i < headerSplit.length; i++) {

                                String colname = headerSplit[i];

                                if (valuesSet.contains(colname)) {

                                    columns.add(colname);
                                    columnIndexes.put(colname, i);

                                }
                            }
                        }

                        headerLine = bean.columns == null ? resultLine
                                : columns.stream()
                                        .collect(
                                                Collectors.joining(IoUtils.separator)
                                        );

                        if (!bean.splitByVariant && !bean.splitByPheno) {

                            String outputPath = bean.outputStem.endsWith(".gz") ? bean.outputStem : bean.outputStem + ".gz";

                            setupWriters(
                                    "generic",
                                    outputPath,
                                    headerComments,
                                    headerLine,
                                    gzWriters,
                                    indexWriters
                            );

                        }

                    } else if (variantIds.isEmpty() && phenoNames.isEmpty()
                            || variantIds.contains(variantId) && phenoNames.isEmpty()
                            || variantIds.isEmpty() && phenoNames.contains(phenoName)
                            || variantIds.contains(variantId) && phenoNames.contains(phenoName)) {

                        String fileKey;
                        if (!bean.splitByPheno && !bean.splitByVariant) {

                            fileKey = "generic";

                        } else if (bean.splitByPheno && !bean.splitByVariant) {

                            fileKey = phenoName;

                        } else if (!bean.splitByPheno && bean.splitByVariant) {

                            fileKey = variantId;

                        } else {

                            fileKey = String.join("_", variantId, phenoName);

                        }

                        IndexedGzWriter outputWriter = gzWriters.get(fileKey);

                        if (outputWriter == null) {

                            String outputPath = String.join(".", bean.outputStem, fileKey, "gz");

                            setupWriters(
                                    fileKey,
                                    outputPath,
                                    headerComments,
                                    headerLine,
                                    gzWriters,
                                    indexWriters
                            );

                        }

                        SimpleFileWriter outputIndexWriter = indexWriters.get(fileKey);

                        String resultLine = gzReader.read(position, compressedLength, uncompressedLength);

                        String newLine;
                        if (bean.columns == null) {

                            newLine = resultLine;

                        } else {

                            String[] resultLineSplit = resultLine.split(IoUtils.separator);

                            newLine = columns.stream()
                                    .mapToInt(
                                            colName -> columnIndexes.get(colName)
                                    )
                                    .mapToObj(
                                            i -> resultLineSplit[i]
                                    )
                                    .collect(
                                            Collectors.joining(IoUtils.separator)
                                    );

                        }

                        IndexedGzCoordinates coordinates = outputWriter.append(newLine);
                        outputIndexWriter.writeLine(
                                variantId,
                                phenoName,
                                Integer.toString(coordinates.compressedLength),
                                Integer.toString(coordinates.uncompressedLength)
                        );

                    }
                }
            }
        } finally {

            gzWriters.values()
                    .forEach(
                            writer -> writer.close()
                    );
            indexWriters.values()
                    .forEach(
                            writer -> writer.close()
                    );

        }
    }

    /**
     * Sets up the writers for the given file key.
     * 
     * @param key the key of the file
     * @param outputPath the path to the output
     * @param headerComments the comments to include in the reader
     * @param headerLine the header line
     * @param gzWriters the map of output gz writers
     * @param indexWriters the map of output index writers
     * 
     * @throws IOException Exception thrown if an I/O error occurs.
     */
    private static void setupWriters(
            String key,
            String outputPath,
            ArrayList<String> headerComments,
            String headerLine,
            HashMap<String, IndexedGzWriter> gzWriters,
            HashMap<String, SimpleFileWriter> indexWriters
    ) throws IOException {

        File outputFile = new File(outputPath);
        IndexedGzWriter outputWriter = new IndexedGzWriter(outputFile);
        gzWriters.put(key, outputWriter);

        File outputIndexFile = IoUtils.getIndexFile(outputFile);
        SimpleFileWriter outputIndexWriter = new SimpleFileWriter(outputIndexFile, true);
        indexWriters.put(key, outputIndexWriter);

        for (String headerComment : headerComments) {

            IndexedGzCoordinates coordinates = outputWriter.append(headerComment);
            outputIndexWriter.writeLine(
                    "Header",
                    "Comment",
                    Integer.toString(coordinates.compressedLength),
                    Integer.toString(coordinates.uncompressedLength)
            );

        }

        IndexedGzCoordinates coordinates = outputWriter.append(headerLine);
        outputIndexWriter.writeLine(
                "Header",
                "Header",
                Integer.toString(coordinates.compressedLength),
                Integer.toString(coordinates.uncompressedLength)
        );
    }

    /**
     * Prints basic help
     */
    private static void printHelp() {

        try (PrintWriter lPrintWriter = new PrintWriter(System.out)) {
            lPrintWriter.print(lineSeparator);
            lPrintWriter.print("==================================" + lineSeparator);
            lPrintWriter.print("              trioGen             " + lineSeparator);
            lPrintWriter.print("               ****               " + lineSeparator);
            lPrintWriter.print("        Results Extraction        " + lineSeparator);
            lPrintWriter.print("==================================" + lineSeparator);
            lPrintWriter.print(lineSeparator
                    + "The extraction command line extracts lines and columns from the results of a linear model." + lineSeparator
                    + lineSeparator
                    + "For documentation and bug report please refer to our code repository https://github.com/mvaudel/trioGen." + lineSeparator
                    + lineSeparator
                    + "----------------------"
                    + lineSeparator
                    + "OPTIONS"
                    + lineSeparator
                    + "----------------------" + lineSeparator
                    + lineSeparator);
            lPrintWriter.print(ExtractOptions.getOptionsAsString());
            lPrintWriter.flush();
        }
    }
}
