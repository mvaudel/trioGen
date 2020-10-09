package no.uib.triogen.io.flat.readers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * Simple wrapper for a flat file reader.
 *
 * @author Marc Vaudel
 */
public class SimpleTextReader implements SimpleFileReader {
    
    /**
     * Lines starting with this character will be ignored.
     */
    private final char COMMENT_CHAR = '#';
    /**
     * Boolean indicating whether comments should be skipped.
     */
    private final boolean skipComments;

    /**
     * The buffered reader.
     */
    private final BufferedReader br;

    /**
     * Constructor.
     *
     * @param file The file to read.
     * @param skipComments Boolean indicating whether comments should be skipped.
     */
    public SimpleTextReader(
            File file,
            boolean skipComments
    ) {
        
        this.skipComments = skipComments;

        try {

            br = new BufferedReader(new FileReader(file));

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    @Override
    public String readLine() {

        try {
            
            String line;
            
            while ((line = br.readLine()) != null) {
                
                if (!skipComments || line.charAt(0) != COMMENT_CHAR) {
                    
                    return line;
                    
                }
            }

            return null;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() {

        try {

            br.close();

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
