package no.uib.triogen.io.variant_id;

import java.io.File;
import java.util.ArrayList;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * Convenience class to read a list of variant ids form a file.
 *
 * @author mvaudel
 */
public class VariantIdFile {

    /**
     * Returns the ids of the variant found in the idsColumn.
     * 
     * @param file The file to read.
     * @param separator The column separator to use.
     * @param idsColumn The name of the column containing the ids.
     * 
     * @return The ids of the variant found in the idsColumn.
     */
    public static String[] parseFromFile(
            File file,
            String separator,
            String idsColumn
    ) {

        ArrayList<String> variantIds = new ArrayList<>();

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(file)) {

            int columnIndex = -1;

            String line = reader.readLine();

            String[] lineSplit = line.split(separator);

            for (int i = 0; i < lineSplit.length; i++) {

                if (lineSplit[i].equals(idsColumn)) {

                    columnIndex = i;
                    break;

                }
            }

            if (columnIndex == -1) {

                try {

                    columnIndex = Integer.parseInt(idsColumn);

                } catch (Exception e) {

                    throw new IllegalArgumentException("Column '" + idsColumn + "' was not found in '" + file + "' using separator '" + separator + "' and cannot be parsed as a number.");

                }

                if (columnIndex < 0) {

                    throw new IllegalArgumentException("Column index '" + idsColumn + "' should be higher or equal to zero.");

                }
                if (columnIndex >= lineSplit.length) {

                    throw new IllegalArgumentException("Column index '" + idsColumn + "' should be 0-based and smaller than the number of columns in the file (" + lineSplit.length + ").");

                }

                variantIds.add(lineSplit[columnIndex]);

            }

            while ((line = reader.readLine()) != null) {

                lineSplit = line.split(separator);
                variantIds.add(lineSplit[columnIndex]);

            }
        }

        return variantIds.stream().toArray(String[]::new);

    }

}
