package no.uib.triogen.model.annotation.ld_link;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import com.mashape.unirest.request.GetRequest;
import java.util.ArrayList;
import static no.uib.triogen.model.annotation.ensembl.EnsemblAPI.getServer;
import no.uib.triogen.model.annotation.ensembl.GeneCoordinates;
import no.uib.triogen.model.annotation.ProxyCoordinates;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * This class queries the LDlink API to retrieve proxies.
 *
 * @author Marc Vaudel
 */
public class LDproxy {

    /**
     * The server url.
     */
    public final static String restUrl = "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy";

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
    public static ArrayList<ProxyCoordinates> getProxy(
            String rsId,
            String pop,
            String r2d,
            String window,
            String token
    ) {

        try {

            String requestUrl = String.join("",
                    restUrl, "?var=", rsId, "&pop=", pop, "&r2_d=", r2d, "&window=", window, "&token=", token
            );
            GetRequest request = Unirest.get(requestUrl);

            HttpResponse<String> answer = request.asString();

            String answerBody = answer.getBody();

            String[] answerBodySplit = answerBody.split("\n");

            ArrayList<ProxyCoordinates> results = new ArrayList<>(answerBodySplit.length - 1);

            String line = answerBodySplit[0];
            String[] lineSplit = line.split("\t");

            if (!lineSplit[0].equals("RS_Number")) {

                throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'RS_Number' was expected in LDproxy query for " + rsId + ".");

            }
            if (!lineSplit[1].equals("Coord")) {

                throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'Coord' was expected in LDproxy query for " + rsId + ".");

            }
            if (!lineSplit[6].equals("R2")) {

                throw new IllegalArgumentException("'" + lineSplit[0] + "' found where 'R2' was expected in LDproxy query for " + rsId + ".");

            }
            if (!lineSplit[7].equals("Correlated_Alleles")) {

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

            return results;

        } catch (Throwable t) {

            throw new RuntimeException(t);

        }
    }
}
