package no.uib.triogen.scripts_marc.ensembl;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.request.GetRequest;
import java.util.ArrayList;
import static no.uib.triogen.model.annotation.EnsemblAPI.getServer;
import no.uib.triogen.model.annotation.VariantCoordinates;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 *
 *
 * @author Marc Vaudel
 */
public class QueryEnsembl {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String rsId = "rs62447179";
        String population = "1000GENOMES:phase_3:GBR";
        double r2Threshold = 0.2;
        
        String ext = String.join("",
                "/ld/homo_sapiens/", rsId, "/", population, "?r2=", Double.toString(r2Threshold), ";attribs=T"
        );
        GetRequest request = Unirest.get(getServer(37) + ext);
        request.header("Content-Type", "application/json");

        try {

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            JSONArray array = jsonResponse.getBody().getArray();
            
            for (int i = 0 ; i < array.length() ; i++) {
                
                JSONObject jsonObject = array.getJSONObject(i);
                
                double r2 = jsonObject.getDouble("r2");
                String proxy = jsonObject.getString("variation");
                int start = jsonObject.getInt("start");
                int end = jsonObject.getInt("end");
            
            int debug = 1;
                
            }

        } catch (Throwable t) {

            t.printStackTrace();

        }

    }

}
