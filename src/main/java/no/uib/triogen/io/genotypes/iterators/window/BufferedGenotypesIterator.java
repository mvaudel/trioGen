package no.uib.triogen.io.genotypes.iterators.window;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.stream.IntStream;
import no.uib.triogen.io.genotypes.GenotypesProvider;
import no.uib.triogen.model.family.ChildToParentMap;
import no.uib.triogen.model.maf.MafEstimator;
import no.uib.triogen.utils.SimpleSemaphore;
import no.uib.triogen.io.genotypes.WindowGenotypesIterator;
import no.uib.triogen.io.genotypes.VariantIterator;
import no.uib.triogen.io.genotypes.iterators.SimpleGenotypeIterator;

/**
 * Buffered iterator for the genotypes of a file.
 *
 * @author Marc Vaudel
 */
public class BufferedGenotypesIterator implements WindowGenotypesIterator {

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
     * Queue of the genotype providers being iterated.
     */
    private final ConcurrentLinkedDeque<GenotypesProvider> currentQueue = new ConcurrentLinkedDeque();
    /**
     * Map of the number of threads per minimal bp needed.
     */
    private final HashMap<String, TreeMap<Integer, Integer>> minBps = new HashMap<>();
    /**
     * Semaphore for the min bp.
     */
    private final SimpleSemaphore minBpSemaphore = new SimpleSemaphore(1);
    /**
     * The number of variants to process in batch.
     */
    private final int nVariants;
    /**
     * Placeholder for a batch of genotypes providers.
     */
    private final ArrayList<GenotypesProvider> batch;
    /**
     * If true, the P0s caches will be filled upon loading and the genotypes
     * provider trimmed.
     */
    private final boolean computeP0;

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
     * @param downstreamLoadingFactor The loading factor to use when trimming
     * the buffer.
     * @param upstreamLoadingFactor The loading factor to use when buffering.
     * @param computeP0 If true, the P0s caches will be filled upon loading and
     * the genotypes provider trimmed.
     */
    public BufferedGenotypesIterator(
            VariantIterator iterator,
            ChildToParentMap childToParentMap,
            int upStreamDistance,
            int downStreamDistance,
            double mafThreshold,
            int nVariants,
            double downstreamLoadingFactor,
            double upstreamLoadingFactor,
            boolean computeP0
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
        this.computeP0 = computeP0;

        init();

    }

    /**
     * Returns the next element. Null if iteration is finished.
     *
     * @return The next element.
     */
    public GenotypesProvider next() {

        if (currentQueue.isEmpty()) {

            bufferSemaphore.acquire();

            bufferSemaphore.release();

            if (!iterator.isFinished()) {

                return next();

            }

            return null;

        }

        GenotypesProvider nextElement = currentQueue.pollFirst();

        String contig = nextElement.getContig();

        buffer(
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

            genotypesProvider.parse(childToParentMap);

            double maf = MafEstimator.getMaf(
                    genotypesProvider,
                    childToParentMap
            );

            if (maf >= mafThreshold) {

                genotypesProvider.setParentP0s(childToParentMap.children, childToParentMap);

                add(genotypesProvider);

                buffer(
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
    private void buffer(
            String contig,
            int bp
    ) {

        if (!iterator.isFinished() && contigList.getLast().equals(contig)) {

            // Load until enough variants are in buffer
            if (bp + upStreamDistance > currentMaxBp.get(contig)) {

                bufferSemaphore.acquire();

                if (!iterator.isFinished() && bp + upStreamDistance > currentMaxBp.get(contig)) {

                    while (!iterator.isFinished() && bp + upstreamLoadingFactor * upStreamDistance >= currentMaxBp.get(contig)) {

                        for (int i = 0; i < 10 * nVariants; i++) {

                            GenotypesProvider genotypesProvider = iterator.next();

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
                                        value -> value.parse(childToParentMap)
                                )
                                .mapToDouble(
                                        value -> MafEstimator.getMaf(
                                                value,
                                                childToParentMap
                                        )
                                )
                                .toArray();

                        if (!batch.isEmpty()) {

                            IntStream.range(0, batch.size())
                                    .parallel()
                                    .filter(
                                            i -> batchMaf[i] >= mafThreshold
                                    )
                                    .mapToObj(
                                            i -> batch.get(i)
                                    )
                                    .forEach(
                                            genotypesProvider -> {

                                                if (computeP0) {

                                                    genotypesProvider.setParentP0s(
                                                            childToParentMap.children,
                                                            childToParentMap
                                                    );

                                                }
                                                
                                                add(genotypesProvider);
                                            
                                            }
                                    );

                            batch.clear();

                        }
                    }

                    // System.out.println("    " + Instant.now() + " - " + bp + " (" + currentMinBp.get(contig) + " - " + currentMaxBp.get(contig) + ", " + buffer.get(contig).size() + " variants in window)");
                }

                bufferSemaphore.release();

            }
        }
    }

    /**
     * Removes data not needed anymore from the buffer.
     *
     * @param contig The contig.
     */
    private void trimBuffer(
            String contig
    ) {

        TreeMap<Integer, Integer> currentMap = minBps.get(contig);

        if (!currentMap.isEmpty()) {

            int minBp = currentMap.firstKey();

            if (minBp - downstreamLoadingFactor * downStreamDistance > currentMinBp.get(contig)) {

                bufferSemaphore.acquire();

                if (!iterator.isFinished() && minBp - downstreamLoadingFactor * downStreamDistance > currentMinBp.get(contig)) {

                    ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>> contigMap = buffer.get(contig);

                    contigMap.keySet().stream()
                            .filter(
                                    tempBp -> tempBp < minBp - downStreamDistance
                            )
                            .forEach(
                                    tempBp -> contigMap.remove(tempBp)
                            );

                    currentMinBp.put(contig, minBp - downStreamDistance);

                    System.gc();

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
    private synchronized void add(
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

            minBpSemaphore.acquire();

            minBps.put(contig, new TreeMap<>());

            minBpSemaphore.release();

        }

        registerMinBp(contig, bp);

        ArrayList<GenotypesProvider> bpMap = contigMap.get(bp);

        if (bpMap == null) {

            bpMap = new ArrayList<>(1);
            contigMap.put(bp, bpMap);

        }

        bpMap.add(genotypesProvider);

    }

    @Override
    public VariantIterator getGenotypesInRange(
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

            bufferSemaphore.acquire();
            bufferSemaphore.release();

            if (!iterator.isFinished() && contigList.getLast().equals(contig)) {

                return getGenotypesInRange(contig, startBp, endBp);

            }
        }

        GenotypesProvider[] genotypes = contigMap.entrySet().stream()
                .filter(
                        entry -> entry.getKey() >= startBp && entry.getKey() <= endBp
                )
                .flatMap(
                        entry -> entry.getValue().stream()
                )
                .toArray(
                        GenotypesProvider[]::new
                );

        return new SimpleGenotypeIterator(genotypes);

    }

    /**
     * Register a lower-bound for the sliding window.
     *
     * @param contig The contig.
     * @param minBp The new lower-bound for the sliding window.
     */
    public void registerMinBp(
            String contig,
            int minBp
    ) {

        boolean removed = false;

        minBpSemaphore.acquire();

        TreeMap<Integer, Integer> currentMap = minBps.get(contig);

        Integer nThreads = currentMap.get(minBp);

        if (nThreads == null) {

            currentMap.put(minBp, 1);

        } else {

            int newValue = nThreads + 1;

            if (newValue != 0) {

                currentMap.put(minBp, newValue);

            } else {

                currentMap.remove(minBp);

                removed = true;

            }
        }

        minBpSemaphore.release();

        if (removed) {

            trimBuffer(contig);

        }
    }

    @Override
    public void releaseMinBp(
            String contig,
            int minBp
    ) {

        boolean removed = false;

        minBpSemaphore.acquire();

        TreeMap<Integer, Integer> currentMap = minBps.get(contig);

        Integer nThreads = currentMap.get(minBp);

        if (nThreads == null) {

            currentMap.put(minBp, -1);

        } else {

            int newValue = nThreads - 1;

            if (newValue == 0) {

                currentMap.remove(minBp);

                removed = true;

                Entry<Integer, Integer> entry = currentMap.firstEntry();

            } else {

                currentMap.put(minBp, newValue);

            }
        }

        minBpSemaphore.release();

        if (removed) {

            trimBuffer(contig);

        }
    }

}
