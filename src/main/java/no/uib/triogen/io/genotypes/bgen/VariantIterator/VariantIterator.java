package no.uib.triogen.io.genotypes.bgen.VariantIterator;

import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.io.genotypes.bgen.reader.BgenFileReader;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Iterator for the variants of a bgen file.
 *
 * @author Marc Vaudel
 */
public class VariantIterator {

    private final BgenIndex bgenIndex;
    private final int start;
    private final int end;

    private final SimpleSemaphore simpleSemaphore = new SimpleSemaphore(1);

    private int currentVariantIndex = 0;

    public VariantIterator(
            BgenIndex bgenIndex,
            int start,
            int end
    ) {

        this.bgenIndex = bgenIndex;
        this.start = start;
        this.end = end;

    }

    public VariantIterator(
            BgenIndex bgenIndex
    ) {

        this(bgenIndex, -1, -1);

    }

    public Integer next() {

        simpleSemaphore.acquire();

        VariantInformation variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

        while (start != -1 && variantInformation.position < start || end != -1 && variantInformation.position > end) {

            currentVariantIndex++;

            if (currentVariantIndex != bgenIndex.variantInformationArray.length) {

                variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

            } else {

                break;

            }
        }

        int index = currentVariantIndex;

        simpleSemaphore.release();

        return index != bgenIndex.variantInformationArray.length ? index : null;

    }

}
