package no.uib.triogen.scripts_marc.cadd;

import io.airlift.compress.zstd.ZstdCompressor;
import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.utils.CompressionUtils;
import no.uib.triogen.utils.TempByteArray;

/**
 * Puts the cadd database in an SQLite database.
 *
 * java -Xmx32G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar
 * no.uib.triogen.scripts_marc.cadd.CaddToSql
 *
 * @author Marc Vaudel
 */
public class CaddToSql {

    public final static int TABLE_SIZE = 65000;
    public final static int BATCH_SIZE = 1000;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File caddFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz");

        File dbFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.sqlite");

        try {

            Connection connection = DriverManager.getConnection("jdbc:sqlite:" + dbFile.getAbsolutePath());
            connection.setAutoCommit(false);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(caddFile, false)) {

                String line = reader.readLine();

                String header = reader
                        .readLine()
                        .substring(1);

                String createStatement = "CREATE TABLE `header` (`id` int, `header` TEXT);";
//        System.out.println(createStatement);
                Statement stmt = connection.createStatement();
                stmt.execute(createStatement);
                connection.commit();

                String insertStatement = "INSERT INTO header (id, header) VALUES (?, ?);";
//        System.out.println(insertStatement);
                PreparedStatement psInsert = connection.prepareStatement(insertStatement);
                psInsert.setInt(1, 0);
                psInsert.setString(2, header);
                psInsert.addBatch();
                psInsert.executeBatch();
                connection.commit();

                String tableColumns = "`id` INTEGER, `chr` TEXT, `bp` INTEGER, `ref` TEXT, `alt` TEXT, `data` TEXT, PRIMARY KEY(id)";
                String headerConcatenated = "id, chr, bp, ref, alt, data";
                String question = "?, ?, ?, ?, ?, ?";

                ArrayList<LineContent> buffer = new ArrayList<>(TABLE_SIZE);

                String lastChromosome = "1";
                int minBp = Integer.MAX_VALUE;
                int maxBp = 0;

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String chromosome = lineSplit[0];
                    int bp = Integer.parseInt(lineSplit[1]);

                    if (bp > maxBp) {

                        maxBp = bp;

                    }
                    if (bp < minBp) {

                        minBp = bp;

                    }

                    byte[] lineBytes = line.getBytes(IoUtils.ENCODING);
                    TempByteArray array = CompressionUtils.zstdCompress(lineBytes);
                    byte[] compressedBytes = Arrays.copyOf(array.array, array.length);

                    String compressedContent = new String(compressedBytes, IoUtils.ENCODING);

                    LineContent lineContent = new LineContent(chromosome, bp, lineSplit[2], lineSplit[3], compressedContent);

                    buffer.add(lineContent);

                    if (!chromosome.equals(lastChromosome) || buffer.size() >= TABLE_SIZE) {

                        writeTable(chromosome, minBp, maxBp, buffer, connection, tableColumns, headerConcatenated, question);

                        lastChromosome = chromosome;
                        buffer.clear();

                    }
                }

                if (!buffer.isEmpty()) {

                    writeTable(lastChromosome, minBp, maxBp, buffer, connection, tableColumns, headerConcatenated, question);

                }
            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }

    private static void writeTable(
            String chromosome,
            int minBp,
            int maxBp,
            ArrayList<LineContent> buffer,
            Connection connection,
            String tableColumns,
            String headerConcatenated,
            String question
    ) throws SQLException {

        String tableName = String.join("_", "table", chromosome, Integer.toString(minBp), Integer.toString(maxBp));

        System.out.println(Instant.now() + " - " + tableName);

        String createStatement = "CREATE TABLE `" + tableName + "` (" + tableColumns + ");";
//        System.out.println(createStatement);
        Statement stmt = connection.createStatement();
        stmt.execute(createStatement);
        connection.commit();

        String insertStatement = "INSERT INTO " + tableName + " (id, " + headerConcatenated + ") VALUES (?, " + question + ");";
//        System.out.println(insertStatement);
        PreparedStatement psInsert = connection.prepareStatement(insertStatement);

        int batchSize = 0;

        for (int i = 0; i < buffer.size(); i++) {

            psInsert.setInt(1, i);

            LineContent lineContent = buffer.get(i);

            psInsert.setString(2, lineContent.chr);
            psInsert.setInt(3, lineContent.bp);
            psInsert.setString(4, lineContent.ref);
            psInsert.setString(5, lineContent.alt);

            psInsert.setString(6, lineContent.compressedData);

            psInsert.addBatch();

            batchSize++;

            if (batchSize >= BATCH_SIZE) {

                psInsert.executeBatch();
                connection.commit();

                psInsert = connection.prepareStatement(insertStatement);

            }

        }

        psInsert.executeBatch();
        connection.commit();

        System.out.println(Instant.now() + " - " + tableName + " Done");

    }

    private static class LineContent {

        public final String chr;
        public final int bp;
        public final String ref;
        public final String alt;
        public final String compressedData;

        public LineContent(
                String chr,
                int bp,
                String ref,
                String alt,
                String compressedData
        ) {
            this.chr = chr;
            this.bp = bp;
            this.ref = ref;
            this.alt = alt;
            this.compressedData = compressedData;
        }
    }

}
