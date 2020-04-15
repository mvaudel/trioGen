package no.uib.triogen.model.annotation;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import com.mashape.unirest.request.GetRequest;
import java.util.ArrayList;
import java.util.HashMap;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * This class allows retrieving information from the Ensembl API.
 *
 * @author Marc Vaudel
 */
public class EnsemblAPI {

    /**
     * The server url.
     */
    public final static String SERVER = "https://rest.ensembl.org";

    /**
     * Returns the coordinates of the features found within the given genomic
     * region as defined as a bp window on a contig.
     *
     * @param contig The contig.
     * @param minBp The minimum of the bp window.
     * @param maxBp The maximum of the bp window.
     *
     * @return The coordinates of the features found within the given genomic
     * region.
     */
    public static ArrayList<GeneCoordinates> getGeneCoordinates(
            String contig,
            int minBp,
            int maxBp
    ) {

        String ext = String.join("",
                "/overlap/region/human/", contig, ":", Integer.toString(minBp), "-", Integer.toString(maxBp), "?feature=gene"
        );
        GetRequest request = Unirest.get(SERVER + ext);
        request.header("Content-Type", "application/json");

        try {

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            JSONArray array = jsonResponse.getBody().getArray();

            ArrayList<GeneCoordinates> result = new ArrayList<>(array.length());

            for (int i = 0; i < array.length(); i++) {

                JSONObject jsonObject = array.getJSONObject(i);

                GeneCoordinates geneCoordinates = new GeneCoordinates(
                        jsonObject.getString("biotype"),
                        jsonObject.has("external_name") ? jsonObject.getString("external_name") : jsonObject.getString("gene_id"),
                        jsonObject.getInt("start"),
                        jsonObject.getInt("end")
                );
                
                result.add(geneCoordinates);

            }

            return result;

        } catch (UnirestException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Returns the Ensembl version as string.
     *
     * @return The Ensembl version as string.
     */
    public static String getEnsemblVersion() {

        String ext = String.join("",
                "/info/data"
        );
        GetRequest request = Unirest.get(SERVER + ext);
        request.header("Content-Type", "application/json");

        try {

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            return jsonResponse.getBody().getObject().get("releases").toString();

        } catch (UnirestException e) {

            throw new RuntimeException(e);

        }

    }

}
