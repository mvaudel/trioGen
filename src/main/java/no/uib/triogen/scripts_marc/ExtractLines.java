package no.uib.triogen.scripts_marc;

import java.io.File;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts the first lines of the given file.
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

        args = new String[]{"C:\\Github\\trioGen\\tmp\\23.vcf.gz", "C:\\Github\\trioGen\\tmp\\23.vcf_firstlines.txt"};

        int start = 0;
        int end = 200;

        try {

            File file = new File(args[0]);
            File destinationFile = new File(args[1]);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(file, false)) {

                try (SimpleFileWriter writer = new SimpleFileWriter(destinationFile, false)) {

                    int count = 0;

                    String line;
                    while ((line = reader.readLine()) != null && ++count <= end) {

                        if (count >= start) {

                            writer.writeLine(line);

                        }
                    }
                }
            }
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }
}
