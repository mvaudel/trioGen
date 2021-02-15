package no.uib.triogen.io.genotypes.bgen.VariantIterator;

import no.uib.triogen.io.genotypes.bgen.index.BgenIndex;
import no.uib.triogen.log.SimpleCliLogger;
import no.uib.triogen.model.genome.VariantInformation;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Iterator for the variants of a bgen file.
 *
 * @author Marc Vaudel
 */
public class VariantIterator {

    private static final int nProgress = 100000;
    private final SimpleCliLogger logger;
    private final String logPrefix;
    private final BgenIndex bgenIndex;
    private final int start;
    private final int end;

    private final SimpleSemaphore simpleSemaphore = new SimpleSemaphore(1);

    private int currentVariantIndex = 0;

    public VariantIterator(
            BgenIndex bgenIndex,
            int start,
            int end,
            SimpleCliLogger logger,
            String logPrefix
    ) {

        this.bgenIndex = bgenIndex;
        this.start = start;
        this.end = end;
        this.logger = logger;
        this.logPrefix = logPrefix;

    }

    public VariantIterator(
            BgenIndex bgenIndex,
            int start,
            int end
    ) {

        this(bgenIndex, start, end, null, null);

    }

    public VariantIterator(
            BgenIndex bgenIndex,
            SimpleCliLogger logger,
            String logPrefix
    ) {

        this(bgenIndex, -1, -1, null, null);

    }

    public VariantIterator(
            BgenIndex bgenIndex
    ) {

        this(bgenIndex, -1, -1, null, null);

    }

    public Integer next() {

        simpleSemaphore.acquire();

        if (logger != null && currentVariantIndex % nProgress == 0) {

            double progress = ((double) (Math.round(10000.0 * currentVariantIndex) / bgenIndex.variantInformationArray.length)) / 100;

            logger.logMessage(logPrefix + "    " + currentVariantIndex + " processed of " + bgenIndex.variantInformationArray.length + " (" + progress + " %)");

        }

        if ((start != -1 || end != -1) && currentVariantIndex < bgenIndex.variantInformationArray.length) {

            VariantInformation variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

            while (start != -1 && variantInformation.position < start || end != -1 && variantInformation.position > end) {

                currentVariantIndex++;

                if (currentVariantIndex != bgenIndex.variantInformationArray.length) {

                    variantInformation = bgenIndex.variantInformationArray[currentVariantIndex];

                } else {

                    break;

                }
            }
        }

        int index = currentVariantIndex;

        currentVariantIndex++;

        simpleSemaphore.release();

        return index < bgenIndex.variantInformationArray.length ? index : null;

    }

}
