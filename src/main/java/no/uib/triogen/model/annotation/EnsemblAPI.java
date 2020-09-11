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
     * The server url for the GRCh38 build.
     */
    public final static String SERVER_GRCh38 = "https://rest.ensembl.org";
    /**
     * The server url for the GRCh37 build.
     */
    public final static String SERVER_GRCh37 = "http://grch37.rest.ensembl.org";

    /**
     * Returns the url of the server to use for the given build.
     *
     * @param buildNumber The number of the build, e.g. 38 for GRCh38.
     *
     * @return The url of the server to use for the given build.
     */
    public static String getServer(
            int buildNumber
    ) {
        if (buildNumber == 37) {

            return SERVER_GRCh37;

        } else if (buildNumber == 38) {

            return SERVER_GRCh38;

        } else {

            throw new UnsupportedOperationException("Build " + buildNumber + " not supported.");

        }
    }

    /**
     * Returns the coordinates of the features found within the given genomic
     * region as defined as a bp window on a contig.
     *
     * @param contig The contig.
     * @param minBp The minimum of the bp window.
     * @param maxBp The maximum of the bp window.
     * @param buildNumber The number of the build, e.g. 38 for GRCh38.
     *
     * @return The coordinates of the features found within the given genomic
     * region.
     */
    public static ArrayList<GeneCoordinates> getGeneCoordinates(
            String contig,
            int minBp,
            int maxBp,
            int buildNumber
    ) {

        String ext = String.join("",
                "/overlap/region/human/", contig, ":", Integer.toString(minBp), "-", Integer.toString(maxBp), "?feature=gene"
        );
        GetRequest request = Unirest.get(getServer(buildNumber) + ext);
        request.header("Content-Type", "application/json");

        try {

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            JSONArray array = jsonResponse.getBody().getArray();

            ArrayList<GeneCoordinates> result = new ArrayList<>(array.length());

            for (int i = 0; i < array.length(); i++) {

                JSONObject jsonObject = array.getJSONObject(i);

                if (jsonObject.has("start") && jsonObject.has("end")) {

                    String bioType = jsonObject.has("biotype") ? jsonObject.getString("biotype") : "Not Available";

                    String geneId = jsonObject.has("external_name")
                            ? jsonObject.getString("external_name")
                            : jsonObject.has("gene_id")
                            ? jsonObject.getString("gene_id")
                            : "Not Available";

                    GeneCoordinates geneCoordinates = new GeneCoordinates(
                            bioType,
                            geneId,
                            jsonObject.getInt("start"),
                            jsonObject.getInt("end")
                    );

                    result.add(geneCoordinates);

                }
            }

            return result;

        } catch (UnirestException e) {

            System.out.println("Faulty Request:\n" + request);
            
            throw new RuntimeException(e);

        }
    }

    /**
     * Returns the Ensembl version as string.
     *
     * @param buildNumber The number of the build, e.g. 38 for GRCh38.
     *
     * @return The Ensembl version as string.
     */
    public static String getEnsemblVersion(
            int buildNumber
    ) {

        String ext = String.join("",
                "/info/data"
        );
        GetRequest request = Unirest.get(getServer(buildNumber) + ext);
        request.header("Content-Type", "application/json");

        try {

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            return jsonResponse.getBody().getObject().get("releases").toString();

        } catch (UnirestException e) {

            throw new RuntimeException(e);

        }

    }

}
