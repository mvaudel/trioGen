package no.uib.triogen.scripts_marc.cadd;

import java.io.File;
import java.time.Instant;
import java.util.ArrayList;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.flat.indexed.IndexedGzCoordinates;
import no.uib.triogen.io.flat.indexed.IndexedGzWriter;

/**
 * Puts the cadd database in an line-indexed gz file.
 *
 * java -Xmx32G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.scripts_marc.cadd.CaddToIndexedGz
 *
 * @author Marc Vaudel
 */
public class CaddToIndexedGz {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File caddFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz");

        File outputFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv.indexed.gz");

        try {

            try (IndexedGzWriter outputWriter = new IndexedGzWriter(outputFile)) {

                SimpleFileWriter index = null;

                try (SimpleFileReader reader = SimpleFileReader.getFileReader(caddFile, false)) {

                    String line = reader.readLine();

                    IndexedGzCoordinates coordinates = outputWriter.append(line + IoUtils.LINE_SEPARATOR);

                    ArrayList<String> indexHeaders = new ArrayList<>(2);

                    String indexHeader = String.join(IoUtils.SEPARATOR,
                            "Header",
                            "Header",
                            "Header",
                            "Comment",
                            Integer.toString(coordinates.compressedLength),
                            Integer.toString(coordinates.uncompressedLength)
                    );
                    indexHeaders.add(indexHeader);

                    line = reader.readLine();
                    coordinates = outputWriter.append(line + IoUtils.LINE_SEPARATOR);

                    indexHeader = String.join(IoUtils.SEPARATOR,
                            "Header",
                            "Header",
                            "Header",
                            "Comment",
                            Integer.toString(coordinates.compressedLength),
                            Integer.toString(coordinates.uncompressedLength)
                    );
                    indexHeaders.add(indexHeader);

                    String lastChromosome = "-1";
                    int lastBp = 0;

                    while ((line = reader.readLine()) != null) {

                        String[] lineSplit = line.split("\t");

                        String chromosome = lineSplit[0];
                        String bpString = lineSplit[1];
                        String ref = lineSplit[2];
                        String alt = lineSplit[3];

                        int bp = Integer.parseInt(bpString);

                        if (bp - lastBp > 100000) {

                            System.out.println(Instant.now() + " -    " + bp);

                            lastBp = bp;

                        }

                        if (!lastChromosome.equals(chromosome)) {

                            System.out.println(Instant.now() + " - Processing Chromosome " + chromosome);

                            File indexFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv." + chromosome + "_index.gz");
                            
                            if (index != null) {
                                
                                index.close();
                                
                            }
                            
                            
                            index = new SimpleFileWriter(indexFile, true);

                            index.writeLine(
                                    "chromosome",
                                    "position",
                                    "ref",
                                    "alt",
                                    "compressedLength",
                                    "uncompressedLength"
                            );

                            for (String headerLine : indexHeaders) {

                                index.writeLine(headerLine);

                            }

                            lastChromosome = chromosome;
                            lastBp = bp;

                        }

                        coordinates = outputWriter.append(line + IoUtils.LINE_SEPARATOR);

                        index.writeLine(
                                chromosome,
                                bpString,
                                ref,
                                alt,
                                Integer.toString(coordinates.compressedLength),
                                Integer.toString(coordinates.uncompressedLength)
                        );
                    }
                }
                
                index.close();
                
            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }

}
