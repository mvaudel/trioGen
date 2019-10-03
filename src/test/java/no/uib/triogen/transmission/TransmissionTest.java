package no.uib.triogen.transmission;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import junit.framework.Assert;
import junit.framework.TestCase;
import no.uib.triogen.cmd.transmission.ExtractTransmission;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This file runs the transmission command on test files and compares the
 * results of the transmission command with ground truth results. Test files are
 * generated by src/R/transmissionTestFiles.R
 *
 * @author Marc Vaudel
 */
public class TransmissionTest extends TestCase {

    /**
     * Runs the command line and checks the output.
     */
    public void testTransmission() {

        HashMap<String, HashMap<String, int[]>> groundTruthMap = getGroundTruth(new File("src/test/resources/transmission/ground_truth.txt"));

        String[] args = new String[]{
            "-g", "src/test/resources/transmission/test_transmission.vcf",
            "-f", "src/test/resources/transmission/test_trio",
            "-o", "src/test/resources/transmission/result.gz"
        };
        ExtractTransmission.main(
                args
        );

        File h1File = new File("src/test/resources/transmission/result.gz_h1.gz");

        if (!h1File.exists()) {

            throw new IllegalArgumentException("Output for h1 not found.");

        }

        HashMap<String, HashMap<String, Integer>> h1Map = getH(h1File);

        File h2File = new File("src/test/resources/transmission/result.gz_h2.gz");

        if (!h2File.exists()) {

            throw new IllegalArgumentException("Output for h2 not found.");

        }

        HashMap<String, HashMap<String, Integer>> h2Map = getH(h2File);

        File h3File = new File("src/test/resources/transmission/result.gz_h3.gz");

        if (!h3File.exists()) {

            throw new IllegalArgumentException("Output for h3 not found.");

        }

        HashMap<String, HashMap<String, Integer>> h3Map = getH(h3File);

        File h4File = new File("src/test/resources/transmission/result.gz_h4.gz");

        if (!h4File.exists()) {

            throw new IllegalArgumentException("Output for h4 not found.");

        }

        HashMap<String, HashMap<String, Integer>> h4Map = getH(h4File);

        for (Entry<String, HashMap<String, int[]>> entry1 : groundTruthMap.entrySet()) {

            String snpId = entry1.getKey();

            for (Entry<String, int[]> entry2 : entry1.getValue().entrySet()) {

                String sampleId = entry2.getKey();
                int[] hs = entry2.getValue();

                if (!h1Map.containsKey(snpId)) {

                    throw new IllegalArgumentException("H1 does not contain " + snpId + ".");

                }
                if (!h1Map.get(snpId).containsKey(sampleId)) {

                    throw new IllegalArgumentException("H1 does not contain sample " + sampleId + " for snp " + snpId + ".");

                }

                int h1Found = h1Map.get(snpId).get(sampleId);
                boolean result = h1Found == hs[0];

                if (!result) {

                    throw new IllegalArgumentException("Incorrect h1 for snp " + snpId + " in sample " + sampleId + ".");

                }

                if (!h2Map.containsKey(snpId)) {

                    throw new IllegalArgumentException("H2 does not contain " + snpId + ".");

                }
                if (!h2Map.get(snpId).containsKey(sampleId)) {

                    throw new IllegalArgumentException("H2 does not contain sample " + sampleId + " for snp " + snpId + ".");

                }

                int h2Found = h2Map.get(snpId).get(sampleId);
                result = h2Found == hs[1];

                if (!result) {

                    throw new IllegalArgumentException("Incorrect h2 for snp " + snpId + " in sample " + sampleId + ".");

                }

                if (!h3Map.containsKey(snpId)) {

                    throw new IllegalArgumentException("H3 does not contain " + snpId + ".");

                }
                if (!h3Map.get(snpId).containsKey(sampleId)) {

                    throw new IllegalArgumentException("H3 does not contain sample " + sampleId + " for snp " + snpId + ".");

                }

                int h3Found = h3Map.get(snpId).get(sampleId);
                result = h3Found == hs[2];

                if (!result) {

                    throw new IllegalArgumentException("Incorrect h3 for snp " + snpId + " in sample " + sampleId + ".");

                }

                if (!h4Map.containsKey(snpId)) {

                    throw new IllegalArgumentException("H4 does not contain " + snpId + ".");

                }
                if (!h4Map.get(snpId).containsKey(sampleId)) {

                    throw new IllegalArgumentException("H4 does not contain sample " + sampleId + " for snp " + snpId + ".");

                }

                int h4Found = h4Map.get(snpId).get(sampleId);
                result = h4Found == hs[3];

                if (!result) {

                    throw new IllegalArgumentException("Incorrect h4 for snp " + snpId + " in sample " + sampleId + ".");

                }

            }
        }
    }

    /**
     * Parses a result file into a map: variant to sample to h.
     *
     * @param hFile the file to parse
     *
     * @return the h
     */
    private HashMap<String, HashMap<String, Integer>> getH(File hFile) {

        HashMap<String, HashMap<String, Integer>> result = new HashMap<>();

        SimpleFileReader reader = SimpleFileReader.getFileReader(hFile);
        String line = reader.readLine();
        String[] header = line.split("\t");

        while ((line = reader.readLine()) != null) {

            String[] lineSplit = line.split("\t");

            String snpId = lineSplit[2];

            HashMap<String, Integer> snpMap = new HashMap<>(header.length - 9);

            for (int i = 8; i < lineSplit.length; i++) {

                String id = header[i];
                int h = Integer.parseInt(lineSplit[i]);
                snpMap.put(id, h);

            }

            result.put(snpId, snpMap);

        }

        return result;

    }

    /**
     * Parses the ground truth results into a map: snpId to sampleId to
     * h1,h2,h3,h4.
     *
     * @param resultFile the result file
     *
     * @return the ground truth results
     */
    private HashMap<String, HashMap<String, int[]>> getGroundTruth(File resultFile) {

        HashMap<String, HashMap<String, int[]>> result = new HashMap<>();

        SimpleFileReader reader = SimpleFileReader.getFileReader(resultFile);

        String line = reader.readLine();

        String[] lineSplit = line.split(" ");
        Assert.assertTrue(lineSplit[0].equals("variant"));
        Assert.assertTrue(lineSplit[1].equals("childId"));
        Assert.assertTrue(lineSplit[10].equals("h1"));
        Assert.assertTrue(lineSplit[11].equals("h2"));
        Assert.assertTrue(lineSplit[12].equals("h3"));
        Assert.assertTrue(lineSplit[13].equals("h4"));

        while ((line = reader.readLine()) != null) {

            lineSplit = line.split(" ");
            String variantId = lineSplit[0];
            String childid = lineSplit[1];
            int h1 = Integer.parseInt(lineSplit[10]);
            int h2 = Integer.parseInt(lineSplit[11]);
            int h3 = Integer.parseInt(lineSplit[12]);
            int h4 = Integer.parseInt(lineSplit[13]);

            HashMap<String, int[]> variantMap = result.get(variantId);

            if (variantMap == null) {

                variantMap = new HashMap<>();
                result.put(variantId, variantMap);

            }

            variantMap.put(childid, new int[]{h1, h2, h3, h4});

        }

        return result;

    }

}
