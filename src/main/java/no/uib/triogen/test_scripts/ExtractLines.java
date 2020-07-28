package no.uib.triogen.test_scripts;

import java.io.File;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts the first lines of the gz files contains in the folder given as argument.
 *
 * @author Marc Vaudel
 */
public class ExtractLines {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        int start = 0;
        int end = 100;

        try {

            File folder = new File(args[0]);

            for (File file : folder.listFiles()) {

                String filePath = file.getAbsolutePath();

                if (filePath.endsWith(".gz")) {

                    File outFile = new File(filePath.substring(0, filePath.length() - 3));

                    try ( SimpleFileReader reader = SimpleFileReader.getFileReader(file)) {

                        try ( SimpleFileWriter writer = new SimpleFileWriter(outFile, false)) {

                            int count = 0;

                            String line;
                            while ((line = reader.readLine()) != null && ++count <= end) {

                                if (count >= start) {

                                    writer.writeLine(line);

                                }
                            }
                        }
                    }
                }
            }
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
