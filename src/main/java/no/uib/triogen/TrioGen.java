package no.uib.triogen;

import java.io.IOException;
import java.io.InputStream;

/**
 *
 * @author Marc Vaudel
 */
public class TrioGen {


    /**
     * Retrieves the version number set in the pom file.
     *
     * @return the version number
     */
    public static String getVersion() {

        java.util.Properties p = new java.util.Properties();

        try {

            InputStream is = (new TrioGen()).getClass().getClassLoader().getResourceAsStream("triogen.properties");
            p.load(is);

        } catch (IOException e) {

            e.printStackTrace();

        }

        return p.getProperty("triogen.version");

//        return "0.5.0-beta";

    }

}
