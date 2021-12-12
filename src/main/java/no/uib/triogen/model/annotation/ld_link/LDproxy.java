package no.uib.triogen.model.annotation.ld_link;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.request.GetRequest;
import java.util.ArrayList;
import no.uib.triogen.model.annotation.ProxyCoordinates;
import no.uib.triogen.utils.Waiter;

/**
 * This class queries the LDlink API to retrieve proxies.
 *
 * @author Marc Vaudel
 */
public class LDproxy {

    /**
     * The server url.
     */
    public final static String REST_URL = "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy";
    /**
     * The default r2/d setting.
     */
    public final static String R2 = "r2";
    /**
     * The default window size.
     */
    public final static String WINDOW = "500000";
    /**
     * The number of attempts.
     */
    public final static int N_ATTEMPTS = 1000;
    /**
     * The delay between attempts.
     */
    public final static long DELAY = 100;
    /**
     * Waiter.
     */
    public final static Waiter WAITER = new Waiter();

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        ArrayList<ProxyCoordinates> debug = getProxy("rs7520742", "CEU&GBR", "r2", "500000", "b823fa70d49d");

        int test = debug.size();

    }

    /**
     * Returns the proxies for the given variant according to LDlink LDproxy.
     *
     * @param rsId The id of the variant.
     * @param pop The population to use.
     * @param r2d A string indicating whether r2 or d should be used.
     * @param window The window to use.
     * @param token The token to use.
     *
     * @return The proxies to use.
     */
    public synchronized static ArrayList<ProxyCoordinates> getProxy(
            String rsId,
            String pop,
            String r2d,
            String window,
            String token
    ) {

        Throwable throwable = null;

        for (int attempt = 0; attempt < N_ATTEMPTS; attempt++) {

            try {

                String requestUrl = String.join("",
                        REST_URL, "?var=", rsId, "&pop=", pop, "&r2_d=", r2d, "&window=", window, "&token=", token
                );
                GetRequest request = Unirest.get(requestUrl);

                HttpResponse<String> answer = request.asString();

                String answerBody = answer.getBody();

                try {

                    ArrayList<ProxyCoordinates> results = new ArrayList<>();

                    if (answerBody.charAt(0) != '{') {

                        String[] answerBodySplit = answerBody.split("\n");

                        String line = answerBodySplit[0];
                        String[] lineSplit = line.split("\t");

                        if (!lineSplit[0].equals("RS_Number")) {

                            System.out.println(answerBody);

                            throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'RS_Number' was expected in LDproxy query for " + rsId + ".");

                        }
                        if (!lineSplit[1].equals("Coord")) {

                            System.out.println(answerBody);

                            throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'Coord' was expected in LDproxy query for " + rsId + ".");

                        }
                        if (!lineSplit[6].equals("R2")) {

                            System.out.println(answerBody);

                            throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'R2' was expected in LDproxy query for " + rsId + ".");

                        }
                        if (!lineSplit[7].equals("Correlated_Alleles")) {

                            System.out.println(answerBody);

                            throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'Correlated_Alleles' was expected in LDproxy query for " + rsId + ".");

                        }

                        for (int lineI = 1; lineI < answerBodySplit.length; lineI++) {

                            line = answerBodySplit[lineI];

                            lineSplit = line.split("\t");

                            String rsid = lineSplit[0];
                            String coord = lineSplit[1];
                            String[] coordSplit = coord.split(":");
                            String chr = coordSplit[0];

                            if (chr.startsWith("chr")) {

                                chr = chr.substring(3);

                            }

                            int bp = Integer.parseInt(coordSplit[1]);
                            double r2 = Double.parseDouble(lineSplit[6]);
                            String correlatedAlleles = lineSplit[7];

                            ProxyCoordinates proxyCoordinates = new ProxyCoordinates(
                                    rsid,
                                    chr,
                                    bp,
                                    bp,
                                    correlatedAlleles,
                                    r2
                            );

                            results.add(proxyCoordinates);

                        }
                    }

                    return results;

                } catch (Throwable t) {

                    System.out.println(answerBody);

                    throwable = t;

                }

            } catch (Throwable t) {

                try {
                    
                    WAITER.delay(DELAY);

                } catch (InterruptedException e) {

                    throw new RuntimeException(e);

                }

                throwable = t;

            }
        }

        throw new RuntimeException(throwable);

    }
}
