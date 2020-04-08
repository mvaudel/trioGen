package no.uib.triogen.io.genotypes.iterators;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedDeque;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.maf.MafEstimator;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * Buffered iterator for the genotypes of a file.
 *
 * @author Marc Vaudel
 */
public class BufferedGenotypesIterator {

    /**
     * The loading factor can be used to reduce the frequency of buffering and
     * cache clean-up. With a loading factor of two, the buffer will be filled
     * twice what is needed.
     */
    private final double downstreamLoadingFactor;
    /**
     * The loading factor can be used to reduce the frequency of buffering and
     * cache clean-up. With a loading factor of two, the buffer will be filled
     * twice what is needed.
     */
    private final double upstreamLoadingFactor;
    /**
     * Highest BP in buffer for each contig.
     */
    private final ConcurrentHashMap<String, Integer> currentMaxBp = new ConcurrentHashMap<>();
    /**
     * Lowest BP in buffer for each contig.
     */
    private final ConcurrentHashMap<String, Integer> currentMinBp = new ConcurrentHashMap<>();
    /**
     * Semaphore for the edition of the buffer.
     */
    private final SimpleSemaphore bufferSemaphore = new SimpleSemaphore(1);
    /**
     * The list of contigs.
     */
    private final ArrayDeque<String> contigList = new ArrayDeque<>();
    /**
     * Contig to bp to snp index to list of array of h.
     */
    private final ConcurrentHashMap<String, ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>>> buffer = new ConcurrentHashMap<>();
    /**
     * The variants iterator.
     */
    private final VariantIterator iterator;
    /**
     * The map of trios.
     */
    private final ChildToParentMap childToParentMap;
    /**
     * The distance in bp to load ahead of the next bp.
     */
    private final int upStreamDistance;
    /**
     * The distance in bp to keep after the next bp.
     */
    private final int downStreamDistance;
    /**
     * The maf threshold. maf is computed in parents and values lower than
     * threshold are not included (inclusive).
     */
    private final double mafThreshold;
    /**
     * The distance in bp to keep after the current bp.
     */
    private final ConcurrentLinkedDeque<GenotypesProvider> currentQueue = new ConcurrentLinkedDeque();
    /**
     * The number of variants to process in batch.
     */
    private final int nVariants;
    /**
     * Placeholder for a batch of genotypes providers.
     */
    private final ArrayList<GenotypesProvider> batch;
    /**
     * Boolean indicating whether the iterator is buffering.
     */
    private boolean buffering = false;

    /**
     * Constructor.
     *
     * @param iterator The variants iterator to use.
     * @param childToParentMap The map of trios.
     * @param upStreamDistance The distance in bp to load ahead of the current
     * bp (inclusive).
     * @param downStreamDistance The distance in bp to keep after the current bp
     * (inclusive).
     * @param mafThreshold The maf threshold. maf is computed in parents and
     * values lower than threshold are not included (inclusive).
     * @param nVariants The number of variants to process in batch.
     * @param downstreamLoadingFactor The loading factor to use when trimming the buffer.
     * @param upstreamLoadingFactor The loading factor to use when buffering.
     */
    public BufferedGenotypesIterator(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            int upStreamDistance,
            int downStreamDistance,
            double mafThreshold,
            int nVariants,
            double downstreamLoadingFactor,
            double upstreamLoadingFactor
    ) {

        this.iterator = iterator;
        this.childToParentMap = childToParentMap;
        this.upStreamDistance = upStreamDistance;
        this.downStreamDistance = downStreamDistance;
        this.mafThreshold = mafThreshold;
        this.nVariants = nVariants;
        this.batch = new ArrayList<>(nVariants);
        this.downstreamLoadingFactor = downstreamLoadingFactor;
        this.upstreamLoadingFactor = upstreamLoadingFactor;
        

        init();

    }

    /**
     * Returns the next element. Null if iteration is finished.
     *
     * @return The next element.
     */
    public GenotypesProvider next() {

        if (currentQueue.isEmpty()) {

            System.out.println("Empty buffer");

            if (buffering) {

                bufferSemaphore.acquire();

                bufferSemaphore.release();

                return next();

            }

            String nextVariant = init();

            if (nextVariant != null) {

                return next();

            }

            return null;

        }

        GenotypesProvider nextElement = currentQueue.pollFirst();

        String contig = nextElement.getContig();

        checkBuffer(
                contig,
                nextElement.getBp()
        );

        while (!contig.equals(contigList.peekFirst())) {

            bufferSemaphore.acquire();

            if (!contig.equals(contigList.peekFirst())) {

                buffer.remove(contigList.pollFirst());

            }

            bufferSemaphore.release();

        }

        return nextElement;

    }

    /**
     * Sets up the iterator.
     *
     * @return Returns the id of the next variant, null if none found.
     */
    private String init() {

        GenotypesProvider genotypesProvider;

        while ((genotypesProvider = iterator.next()) != null) {

            genotypesProvider.parse();

            double maf = MafEstimator.getMaf(
                    genotypesProvider,
                    childToParentMap
            );

            if (maf >= mafThreshold) {

                add(genotypesProvider);

                checkBuffer(
                        genotypesProvider.getContig(),
                        genotypesProvider.getBp()
                );

                return genotypesProvider.getVariantID();

            }
        }

        return null;

    }

    /**
     * Checks that the buffer contains enough data for the given contig at the
     * given bp.
     *
     * @param contig The contig to buffer.
     * @param bp The current bp.
     */
    private void checkBuffer(
            String contig,
            int bp
    ) {

        if (contigList.getLast().equals(contig)) {

            // Load until enough variants are in buffer
            if (bp + upStreamDistance > currentMaxBp.get(contig)) {

                bufferSemaphore.acquire();

                if (bp + upStreamDistance > currentMaxBp.get(contig)) {

                    buffering = true;

                    while (bp + upstreamLoadingFactor * upStreamDistance >= currentMaxBp.get(contig)) {

                        GenotypesProvider genotypesProvider;

                        for (int i = 0; i < nVariants; i++) {

                            genotypesProvider = iterator.next();

                            if (genotypesProvider == null) {

                                break;

                            }

                            batch.add(genotypesProvider);

                        }

                        if (batch.isEmpty()) {

                            return;

                        }

                        double[] batchMaf = batch.parallelStream()
                                .peek(
                                        value -> value.parse()
                                )
                                .mapToDouble(
                                        value -> MafEstimator.getMaf(
                                                value,
                                                childToParentMap
                                        )
                                )
                                .toArray();

                        for (int i = 0; i < batch.size(); i++) {

                            if (batchMaf[i] >= mafThreshold) {

                                add(batch.get(i));

                            }
                        }

                        batch.clear();

                    }

                    buffering = false;

                    System.out.println("New window: " + bp + " (" + currentMinBp.get(contig) + " - " + currentMaxBp.get(contig) + ", " + buffer.get(contig).size() + " variants in buffer)");

                }

                bufferSemaphore.release();

            }

            // Remove variants outside range
            if (bp - downstreamLoadingFactor * downStreamDistance > currentMinBp.get(contig)) {

                bufferSemaphore.acquire();

                if (bp - downstreamLoadingFactor * downStreamDistance > currentMinBp.get(contig)) {

                    ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>> contigMap = buffer.get(contig);

                    contigMap.keySet().stream()
                            .filter(
                                    tempBp -> tempBp < bp - downStreamDistance
                            )
                            .forEach(
                                    tempBp -> contigMap.remove(tempBp)
                            );

                    currentMinBp.put(contig, bp - downStreamDistance);

                    System.gc();

                    System.out.println("New window: " + bp + " (" + currentMinBp.get(contig) + " - " + currentMaxBp.get(contig) + ", " + buffer.get(contig).size() + " variants in buffer)");

                }

                bufferSemaphore.release();

            }
        }
    }

    /**
     * Adds a variant to the buffer.
     *
     * @param genotypesProvider The genotypes provider for this variant.
     */
    private void add(
            GenotypesProvider genotypesProvider
    ) {

        String contig = genotypesProvider.getContig();
        int bp = genotypesProvider.getBp();

        currentQueue.add(genotypesProvider);

        currentMaxBp.put(contig, bp);

        ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>> contigMap = buffer.get(contig);

        if (contigMap == null) {

            contigMap = new ConcurrentHashMap<>();
            buffer.put(contig, contigMap);
            contigList.add(contig);

            currentMinBp.put(contig, 0);

        }

        ArrayList<GenotypesProvider> bpMap = contigMap.get(bp);

        if (bpMap == null) {

            bpMap = new ArrayList<>(1);
            contigMap.put(bp, bpMap);

        }

        bpMap.add(genotypesProvider);

    }

    /**
     * Returns an array of the genotypes of the given contig in the given bp
     * range.
     *
     * @param contig The contig.
     * @param startBp The bp range start (inclusive).
     * @param endBp The bp range end (inclusive).
     *
     * @return An array of the genotypes.
     */
    public GenotypesProvider[] getGenotypesInRange(
            String contig,
            int startBp,
            int endBp
    ) {

        ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>> contigMap = buffer.get(contig);

        if (contigMap == null) {

            throw new IllegalArgumentException("Contig " + contig + " not in buffer.");

        }

        if (startBp < currentMinBp.get(contig)) {

            throw new IllegalArgumentException("Sliding window start (" + startBp + ") out of range (" + currentMinBp.get(contig) + " - " + currentMaxBp.get(contig) + ").");

        }

        if (endBp > currentMaxBp.get(contig)) {

            throw new IllegalArgumentException("Sliding window end (" + endBp + ") out of range (" + currentMinBp.get(contig) + " - " + currentMaxBp.get(contig) + ").");

        }

        return contigMap.entrySet().stream()
                .filter(
                        entry -> entry.getKey() >= startBp && entry.getKey() <= endBp
                )
                .flatMap(
                        entry -> entry.getValue().stream()
                )
                .toArray(
                        GenotypesProvider[]::new
                );
    }

}
