package no.uib.triogen.processing.ld;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Cache for p0 values.
 *
 * @author Marc Vaudel
 */
public class P0Cache {

    private final int[] checkOutPosition;

    private final TreeMap<Integer, ArrayList<String>> positionMap = new TreeMap<>();

    private final ConcurrentHashMap<String, float[][]> pHomozygous = new ConcurrentHashMap<>();

    public P0Cache(
            int nThreads
    ) {

        this.checkOutPosition = new int[nThreads];

    }

    public float[][] getPHomozygous(
            String id
    ) {

        return pHomozygous.get(id);

    }

    public void release(int thread, int pos) {

        checkOutPosition[thread] = pos;

        for (int threadPosition : checkOutPosition) {

            if (threadPosition < pos) {

                return;

            }
        }

        for (Entry<Integer, ArrayList<String>> entry : positionMap.entrySet()) {

            int position = entry.getKey();

            if (position > pos) {

                return;

            }

            for (String id : entry.getValue()) {

                pHomozygous.remove(id);

            }
        }
    }

    public void register(
            BgenVariantData variantData,
            ChildToParentMap childToParentMap
    ) {

        VariantInformation variantInformation = variantData.getVariantInformation();

        float[][] variantPHomozygous = new float[variantInformation.alleles.length][2 * childToParentMap.children.length];

        for (int alleleI = 0; alleleI <= variantInformation.alleles.length; alleleI++) {

            for (int i = 0; i < childToParentMap.children.length; i++) {

                String childId = childToParentMap.children[i];

                String motherId = childToParentMap.getMother(childId);

                if (variantData.contains(motherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(motherId); z++) {

                        homozygous *= variantData.getProbability(motherId, z, alleleI);

                    }

                    variantPHomozygous[alleleI][i] = homozygous;

                } else {

                    variantPHomozygous[alleleI][i] = Float.NaN;

                }

                String fatherId = childToParentMap.getFather(childId);

                if (variantData.contains(motherId)) {

                    float homozygous = 1.0f;

                    for (int z = 0; z < variantData.getPloidy(fatherId); z++) {

                        homozygous *= variantData.getProbability(fatherId, z, alleleI);

                    }

                    variantPHomozygous[alleleI][i + childToParentMap.children.length] = homozygous;

                } else {

                    variantPHomozygous[alleleI][i + childToParentMap.children.length] = Float.NaN;

                }
            }
        }

        pHomozygous.put(variantInformation.id, variantPHomozygous);

        ArrayList<String> idsAtPosition = positionMap.get(variantInformation.position);

        if (idsAtPosition == null) {

            idsAtPosition = new ArrayList<>(2);
            positionMap.put(variantInformation.position, idsAtPosition);

        }

        idsAtPosition.add(variantInformation.id);

    }
}
