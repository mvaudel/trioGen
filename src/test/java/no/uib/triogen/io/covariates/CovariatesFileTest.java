package no.uib.triogen.io.covariates;

import java.util.HashMap;
import java.util.TreeSet;
import junit.framework.Assert;
import junit.framework.TestCase;

/**
 *
 *
 * @author Marc Vaudel
 */
public class CovariatesFileTest extends TestCase {

    public void testParsing() {

        String fileContent = "{\"z_pregnancy_duration\":{\"indivudual\":[\"child\"],\"phenoName\":[\"pregnancy_duration\"],\"phenoLabel\":[\"Pregnancy duration\"],\"phenoUnit\":[\"Day\"],\"formula\":[\"pregnancy_duration ~ 1\"],\"sigmaFormula\":[\"~ 1\"],\"family\":[\"BCT\"],\"controlVariables\":[\"sex_number\"],\"controlVariablesLabels\":[\"Sex\"],\"controlVariablesUnits\":[\"1: male, 2: female\"],\"controlVariablesType\":[\"discrete\"]},\"z_umbilical_cord_length\":{\"indivudual\":[\"child\"],\"phenoName\":[\"umbilical_cord_length\"],\"phenoLabel\":[\"Umbilical cord length\"],\"phenoUnit\":[\"cm\"],\"formula\":[\"umbilical_cord_length ~ fp(pregnancy_duration)\"],\"sigmaFormula\":[\" ~ fp(pregnancy_duration)\"],\"family\":[\"BCT\"],\"controlVariables\":[\"sex_number\",\"pregnancy_duration\"],\"controlVariablesLabels\":[\"Sex\",\"Pregnancy duration\"],\"controlVariablesUnits\":[\"1: male, 2: female\",\"Day\"],\"controlVariablesType\":[\"discrete\",\"continuous\"]},\"z_placenta_weight\":{\"indivudual\":[\"child\"],\"phenoName\":[\"placenta_weight\"],\"phenoLabel\":[\"Placenta weight\"],\"phenoUnit\":[\"g\"],\"formula\":[\"placenta_weight ~ fp(pregnancy_duration)\"],\"sigmaFormula\":[\" ~ fp(pregnancy_duration)\"],\"family\":[\"BCT\"],\"controlVariables\":[\"sex_number\",\"pregnancy_duration\"],\"controlVariablesLabels\":[\"Sex\",\"Pregnancy duration\"],\"controlVariablesUnits\":[\"1: male, 2: female\",\"Day\"],\"controlVariablesType\":[\"discrete\",\"continuous\"]},\"z_no_control_variable\":{\"indivudual\":[\"child\"],\"phenoName\":[\"no_control_variable\"],\"phenoLabel\":[\"Mother diet alcohol intake\"],\"phenoUnit\":[\"Uknown unit\"],\"formula\":[\"no_control_variable ~ none\"],\"sigmaFormula\":[\" ~ none\"],\"family\":[\"NO\"],\"controlVariables\":{},\"controlVariablesLabels\":{},\"controlVariablesUnits\":{},\"controlVariablesType\":{}}}";

        HashMap<String, TreeSet<String>> covariatesMap = SpecificCovariatesFile.praseCovariates(fileContent);

        Assert.assertTrue(covariatesMap.size() == 4);

        Assert.assertTrue(covariatesMap.containsKey("z_pregnancy_duration"));
        TreeSet<String> covariates = covariatesMap.get("z_pregnancy_duration");
        Assert.assertTrue(covariates.size() == 1);
        Assert.assertTrue(covariates.contains("sex_number"));

        Assert.assertTrue(covariatesMap.containsKey("z_umbilical_cord_length"));
        covariates = covariatesMap.get("z_umbilical_cord_length");
        Assert.assertTrue(covariates.size() == 2);
        Assert.assertTrue(covariates.contains("sex_number"));
        Assert.assertTrue(covariates.contains("pregnancy_duration"));

        Assert.assertTrue(covariatesMap.containsKey("z_placenta_weight"));
        covariates = covariatesMap.get("z_placenta_weight");
        Assert.assertTrue(covariates.size() == 2);
        Assert.assertTrue(covariates.contains("sex_number"));
        Assert.assertTrue(covariates.contains("pregnancy_duration"));

        Assert.assertTrue(covariatesMap.containsKey("z_no_control_variable"));
        covariates = covariatesMap.get("z_no_control_variable");
        Assert.assertTrue(covariates.isEmpty());

        Assert.assertTrue(!covariatesMap.containsKey("dummy_key"));

    }
}
