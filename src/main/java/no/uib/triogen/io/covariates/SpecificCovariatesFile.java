package no.uib.triogen.io.covariates;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * Convenience class used to parse covariates specification files.
 *
 * @author Marc Vaudel
 */
public class SpecificCovariatesFile {

    /**
     * Parses the content of a covariates specification file. Note: IO
     * exceptions encountered while reading the file are thrown as runtime
     * exceptions.
     *
     * @param file The file to read.
     *
     * @return The covariates to use for every phenotype.
     */
    public static HashMap<String, TreeSet<String>> praseCovariates(
            File file
    ) {

        try {

            String fileContent = Files.lines(file.toPath())
                    .collect(
                            Collectors.joining()
                    );

            return praseCovariates(fileContent);

        } catch (IOException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Parses the content of a covariates specification file.
     *
     * @param fileContent The content of the file as a String.
     *
     * @return The covariates to use for every phenotype.
     */
    public static HashMap<String, TreeSet<String>> praseCovariates(
            String fileContent
    ) {

        JSONObject jsonObject = new JSONObject(fileContent);

        Set<String> keys = jsonObject.keySet();

        HashMap<String, TreeSet<String>> result = new HashMap<>(keys.size());

        for (String key : keys) {

            JSONObject phenoObject = jsonObject.getJSONObject(key);

            TreeSet<String> covariates = new TreeSet<>();

            result.put(key, covariates);

            Object controlVariablesObject = phenoObject.get("controlVariables");

            if (controlVariablesObject instanceof JSONArray) {

                JSONArray array = (JSONArray) controlVariablesObject;

                for (int i = 0; i < array.length(); i++) {

                    covariates.add(array.getString(i));

                }
            }
        }

        return result;

    }

}
