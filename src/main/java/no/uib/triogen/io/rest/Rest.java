package no.uib.triogen.io.rest;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.request.GetRequest;
import no.uib.triogen.utils.SimpleTimeSemaphore;

/**
 * Functions to interact with a rest API.
 *
 * @author Marc Vaudel
 */
public class Rest {
    
    /**
     * Simple time semaphore to increase delay between queries.
     */
    private static SimpleTimeSemaphore simpleTimeSemaphore = new SimpleTimeSemaphore();

    /**
     * Returns the response as json.
     * 
     * @param request The request.
     * 
     * @return The response as json.
     */
    public static JsonNode getJson(GetRequest request) {

        return getResponse(request).getBody();

    }

    /**
     * Returns the response.
     * 
     * @param request The request.
     * 
     * @return The response.
     */
    public static HttpResponse<JsonNode> getResponse(
            GetRequest request
    ) {
        
        return getResponse(request, 0);
        
    }

    /**
     * Returns the response.
     * 
     * @param request The request.
     * @param attempts The number of prior attempts.
     * 
     * @return The response.
     */
    private static HttpResponse<JsonNode> getResponse(
            GetRequest request, 
            int attempts
    ) {

        try {

            simpleTimeSemaphore.delay();

            HttpResponse<JsonNode> jsonResponse = request.asJson();

            if (jsonResponse.getStatus() == 429) {

                simpleTimeSemaphore.increaseDelay();

                return getResponse(request);

            } else if (jsonResponse.getStatus() != 200) {

                throw new RuntimeException("HTTP error: "
                        + jsonResponse.getStatus()
                        + " - "
                        + jsonResponse.getStatusText());
            }

            return jsonResponse;

        } catch (Exception e) {
            
            if (attempts++ < 100) {
                
                return getResponse(request, attempts + 1);
                
            }

            try {

                System.out.println(request.asString().getBody());

            } catch (Exception e2) {
                e2.printStackTrace();
            }

            throw new RuntimeException(e);
        }

    }
}
