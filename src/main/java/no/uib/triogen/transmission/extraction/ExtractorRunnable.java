package no.uib.triogen.transmission.extraction;

import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import no.uib.triogen.io.flat.SimpleFileWriter;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.model.family.ChildToParentMap;

/**
 * Runnable for the extraction of transmitted alleles.
 *
 * @author Marc Vaudel
 */
public class ExtractorRunnable implements Runnable {

    /**
     * The variants iterator.
     */
    private final VariantIterator iterator;
    /**
     * The h1 writer.
     */
    private final SimpleFileWriter h1Writer;
    /**
     * The h2 writer.
     */
    private final SimpleFileWriter h2Writer;
    /**
     * The h3 writer.
     */
    private final SimpleFileWriter h3Writer;
    /**
     * The h4 writer.
     */
    private final SimpleFileWriter h4Writer;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * Boolean indicating whether the runnable has been canceled.
     */
    private static boolean canceled = false;

    /**
     * Constructor.
     *
     * @param iterator the variants iterator
     * @param childToParentMap the child to parent map
     * @param h1Writer the h1 writer
     * @param h2Writer the h2 writer
     * @param h3Writer the h3 writer
     * @param h4Writer the h4 writer
     */
    public ExtractorRunnable(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            SimpleFileWriter h1Writer,
            SimpleFileWriter h2Writer,
            SimpleFileWriter h3Writer,
            SimpleFileWriter h4Writer
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.h1Writer = h1Writer;
        this.h2Writer = h2Writer;
        this.h3Writer = h3Writer;
        this.h4Writer = h4Writer;

    }

    @Override
    public void run() {

        try {

            GenotypesProvider genotypesProvider;
            while ((genotypesProvider = iterator.next()) != null && !canceled) {

                genotypesProvider.parse();

                String[] genotypes = extractHs(genotypesProvider);

                h1Writer.writeLine(
                        String.join(
                                "\t",
                                genotypesProvider.getVariantDescription(),
                                genotypes[0]
                        )
                );
                h2Writer.writeLine(
                        String.join(
                                "\t",
                                genotypesProvider.getVariantDescription(),
                                genotypes[1]
                        )
                );
                h3Writer.writeLine(
                        String.join(
                                "\t",
                                genotypesProvider.getVariantDescription(),
                                genotypes[2]
                        )
                );
                h4Writer.writeLine(
                        String.join(
                                "\t",
                                genotypesProvider.getVariantDescription(),
                                genotypes[3]
                        )
                );

            }
        } catch (Throwable t) {

            canceled = true;
            t.printStackTrace();

        }
    }

    /**
     * Extracts the Hs from the given genotypes.
     *
     * @param genotypesProvider a genotypes provider
     *
     * @return an array of the different hs to export
     */
    private String[] extractHs(
            GenotypesProvider genotypesProvider
    ) {

        HashMap<String, int[]> hMap = childToParentMap.children.stream()
                .parallel()
                .collect(
                        Collectors.toMap(
                                childId -> childId,
                                childId -> getH(
                                        genotypesProvider,
                                        childId
                                ),
                                (a, b) -> a,
                                HashMap::new
                        )
                );

        return IntStream.range(0, 4)
                .parallel()
                .mapToObj(
                        i -> aggregateH(hMap, i)
                )
                .toArray(
                        String[]::new
                );

    }

    /**
     * Aggregates an array of h in a string.
     *
     * @param hs the h matrix
     * @param i the index of h to aggregate
     *
     * @return a tab separated string of the hs of all kids
     */
    private String aggregateH(
            HashMap<String, int[]> hMap,
            int i
    ) {

        StringBuilder sb = new StringBuilder(2 * hMap.size() - 1);

        Iterator<String> childIt = childToParentMap.children.iterator();

        sb.append(hMap.get(childIt.next())[i]);

        while (childIt.hasNext()) {

            sb.append('\t')
                    .append(hMap.get(childIt.next())[i]);

        }

        return sb.toString();

    }

    /**
     * Returns an array containing h1, h2, h3, and h4 for a given child using the given genotypes provider.
     *
     * @param vcfLine the vcf line
     * @param childId the child id
     *
     * @return an array containing h1, h2, h3, and h4
     */
    private int[] getH(
            GenotypesProvider genotypesProvider,
            String childId
    ) {

        String motherId = childToParentMap.getMother(childId);
        String fatherId = childToParentMap.getFather(childId);

        int genotypeKid = genotypesProvider.getGenotype(childId);
        int genotypeMother = genotypesProvider.getGenotype(motherId);
        int genotypeFather = genotypesProvider.getGenotype(fatherId);

        int nAltMother = genotypeMother >= 2 ? genotypeMother - 1 : genotypeMother;
        int nAltFather = genotypeFather >= 2 ? genotypeFather - 1 : genotypeFather;

        int h3 = genotypeKid == 0 || genotypeKid == 2 ? 0 : 1;
        int h1 = genotypeKid == 0 || genotypeKid == 1 ? 0 : 1;
        int h2 = nAltMother - h1;
        int h4 = nAltFather - h3;

        return new int[]{h1, h2, h3, h4};

    }
}
