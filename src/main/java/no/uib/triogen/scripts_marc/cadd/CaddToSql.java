package no.uib.triogen.scripts_marc.cadd;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * Puts the cadd database in an SQLite database.
 * 
 * java -Xmx32G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.scripts_marc.cadd.CaddToSql
 *
 * @author Marc Vaudel
 */
public class CaddToSql {

    public final static int TABLE_SIZE = 65000;

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

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(caddFile)) {

                String line = reader.readLine();

                String[] header = reader
                        .readLine()
                        .substring(1)
                        .split("\t");

                String[] headerReformatted = Arrays.stream(header)
                        .map(
                                content -> content
                                        .replace('-', '_')
                                        .replace(' ', '_')
                                        .replace(',', '_')
                        )
                        .toArray(
                                String[]::new
                        );

                String headerConcatenated = String.join(", ", headerReformatted);

                String question = Arrays.stream(header)
                        .map(
                                content -> "?"
                        )
                        .collect(
                                Collectors.joining(", ")
                        );

                StringBuilder stringBuilder = new StringBuilder("`id` INTEGER");

                for (String colName : headerReformatted) {

                    stringBuilder.append(", `")
                            .append(colName)
                            .append("` TEXT");

                }

                stringBuilder.append(", PRIMARY KEY(id)");

                String tableColumns = stringBuilder.toString();

                HashMap<Integer, String[]> buffer = new HashMap<>(TABLE_SIZE);

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

                    int id = line.hashCode();

                    if (!buffer.containsKey(id)) {

                        buffer.put(id, lineSplit);

                    } else {

                        throw new IllegalArgumentException("Duplicate entry: '" + id + "'.");

                    }

                    if (!chromosome.equals(lastChromosome) || buffer.size() >= TABLE_SIZE) {

                        writeTable(chromosome, minBp, maxBp, buffer, connection, tableColumns, headerConcatenated, question);

                        lastChromosome = chromosome;
                        buffer.clear();

                    }
                }

                writeTable(lastChromosome, minBp, maxBp, buffer, connection, tableColumns, headerConcatenated, question);

            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }

    private static void writeTable(
            String chromosome,
            int minBp,
            int maxBp,
            HashMap<Integer, String[]> buffer,
            Connection connection,
            String tableColumns,
            String headerConcatenated,
            String question
    ) throws SQLException {

        String tableName = String.join("_", chromosome, Integer.toString(minBp), Integer.toString(maxBp));

        System.out.println(Instant.now() + " - " + tableName);

        String sql = "CREATE TABLE `" + tableName + "` (" + tableColumns + ");";
        Statement stmt = connection.createStatement();
        stmt.execute(sql);

        String insertStatement = "INSERT INTO " + tableName + " (id, " + headerConcatenated + ") VALUES (?, " + question + ");";
        PreparedStatement psInsert = connection.prepareStatement(insertStatement);

        for (Entry<Integer, String[]> entry : buffer.entrySet()) {

            psInsert.setInt(1, entry.getKey());

            for (int i = 0; i < entry.getValue().length; i++) {

                psInsert.setString(i + 2, entry.getValue()[i]);

            }

            psInsert.addBatch();

        }

        psInsert.executeBatch();
        connection.commit();

    }

}
