package no.uib.triogen.io.flat.readers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.zip.GZIPInputStream;
import no.uib.triogen.io.flat.SimpleFileReader;
import static no.uib.triogen.io.IoUtils.ENCODING;

/**
 * Simple wrapper for a gz file reader.
 *
 * @author Marc Vaudel
 */
public class SimpleGzReader implements SimpleFileReader {
    
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
    public SimpleGzReader(
            File file,
            boolean skipComments
    ) {
        
        this.skipComments = skipComments;

        try {

            InputStream fileStream = new FileInputStream(file);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, ENCODING);

            br = new BufferedReader(decoder);

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
