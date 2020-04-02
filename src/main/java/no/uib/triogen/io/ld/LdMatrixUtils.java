package no.uib.triogen.io.ld;

import java.io.UnsupportedEncodingException;
import static no.uib.triogen.io.IoUtils.ENCODING;

/**
 * Utils for the reading and writing of ld matrices.
 *
 * @author Marc Vaudel
 */
public class LdMatrixUtils {

    /**
     * The file extension.
     */
    public static final String EXTENSION = ".tld";
    /**
     * The magic number of to use to identify the supported files.
     */
    public static final byte[] MAGIC_NUMBER = getMagicNumber();

    /**
     * Returns the magic number.
     *
     * @return The magic number.
     */
    public static byte[] getMagicNumber() {

        try {

            String magicName = "Triogen.ldMatrix.1";
            return magicName.getBytes(ENCODING);

        } catch (UnsupportedEncodingException e) {

            throw new RuntimeException(e);

        }
    }

}
