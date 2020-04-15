package no.uib.triogen.test_scripts;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.request.GetRequest;
import no.uib.triogen.model.annotation.EnsemblAPI;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * Sand box for the Ensembl API code.
 *
 * @author Marc Vaudel
 */
public class EnsemblAPITest {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {

            String version = EnsemblAPI.getEnsemblVersion();

        } catch (Throwable t) {
            t.printStackTrace();
        }
    }

}
