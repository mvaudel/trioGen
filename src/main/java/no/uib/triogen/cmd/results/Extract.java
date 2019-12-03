package no.uib.triogen.cmd.results;

import java.io.File;
import no.uib.triogen.cmd.association.*;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.lineSeparator;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.vcf.custom.CustomVcfIterator;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.geno.Model;
import no.uib.triogen.model.geno.VariantList;
import no.uib.triogen.processing.association.linear_model.LinearModelComputer;
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
     */
    private static void run(
            ExtractOptionsBean bean
    ) {

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(bean.inputFile)) {

            String line = reader.readLine();
            String[] headerSplit = line.split(IoUtils.separator);

            ArrayList<String> values = new ArrayList<>(headerSplit.length);
            ArrayList<String> categories = new ArrayList<>(1);
            HashMap<String, Integer> columnIndexes = new HashMap<>();

            if (bean.valueColumns != null || bean.categoryColumns != null) {

                HashSet<String> valuesSet = bean.valueColumns == null ? new HashSet<>(0)
                        : Arrays.stream(bean.valueColumns)
                                .collect(
                                        Collectors.toCollection(HashSet::new)
                                );

                HashSet<String> categoriesSet = bean.categoryColumns == null ? new HashSet<>(0)
                        : Arrays.stream(bean.categoryColumns)
                                .collect(
                                        Collectors.toCollection(HashSet::new)
                                );

                HashSet<String> headerSet = Arrays.stream(headerSplit)
                        .collect(
                                Collectors.toCollection(HashSet::new)
                        );

                String missingColumns = Stream.concat(valuesSet.stream(), categoriesSet.stream())
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

                        values.add(colname);
                        columnIndexes.put(colname, i);

                    }
                    if (categoriesSet.contains(colname)) {

                        categories.add(colname);
                        columnIndexes.put(colname, i);

                    }
                }
            }

            String headerLine = bean.valueColumns == null ? line
                    : values.stream()
                            .collect(
                                    Collectors.joining(IoUtils.separator)
                            );

            HashMap<String, SimpleFileWriter> fileWriters = new HashMap<>();

            if (bean.categoryColumns == null) {

                String outputPath = bean.outputStem.endsWith(".gz") ? bean.outputStem : bean.outputStem + ".gz";
                SimpleFileWriter writer = new SimpleFileWriter(new File(outputPath), true);
                writer.writeLine(headerLine);
                fileWriters.put("dummy", writer);

            }

            while ((line = reader.readLine()) != null) {

                String[] lineSplit = line.split(IoUtils.separator);

                String phenoKey = bean.categoryColumns == null
                        ? "dummy"
                        : categories.stream()
                                .mapToInt(
                                        colName -> columnIndexes.get(colName)
                                )
                                .mapToObj(
                                        i -> lineSplit[i]
                                )
                                .collect(
                                        Collectors.joining("_")
                                );

                SimpleFileWriter writer = fileWriters.get(phenoKey);

                if (writer == null) {

                    String outputPath = String.join(".", bean.outputStem, phenoKey, "gz");
                    writer = new SimpleFileWriter(new File(outputPath), true);
                    writer.writeLine(headerLine);
                    fileWriters.put(phenoKey, writer);

                }

                String newLine = bean.valueColumns == null ? line
                        : values.stream()
                                .mapToInt(
                                        colName -> columnIndexes.get(colName)
                                )
                                .mapToObj(
                                        i -> lineSplit[i]
                                )
                                .collect(
                                        Collectors.joining(IoUtils.separator)
                                );

                writer.writeLine(newLine);

            }
            
            fileWriters.values()
                    .forEach(
                            writer -> writer.close()
                    );
            
        }
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
            lPrintWriter.print("   Transmitted Allele Extraction  " + lineSeparator);
            lPrintWriter.print("==================================" + lineSeparator);
            lPrintWriter.print(lineSeparator
                    + "The extraction command line iterates a vcf file and extracts transmitted alleles of children." + lineSeparator
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
