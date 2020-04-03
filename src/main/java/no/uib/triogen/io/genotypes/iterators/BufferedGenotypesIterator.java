package no.uib.triogen.io.genotypes.iterators;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedDeque;
import no.uib.triogen.io.genotypes.GenotypesProvider;
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
    public static final double LOADING_FACTOR = 2.0;
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
     * The distance in bp to load ahead of the next bp.
     */
    private final int upStreamDistance;
    /**
     * The distance in bp to keep after the next bp.
     */
    private final int downStreamDistance;
    /**
     * The distance in bp to keep after the current bp.
     */
    private final ConcurrentLinkedDeque<GenotypesProvider> currentQueue = new ConcurrentLinkedDeque();

    /**
     * Constructor.
     *
     * @param iterator The variants iterator to use.
     * @param upStreamDistance The distance in bp to load ahead of the current
     * bp (inclusive).
     * @param downStreamDistance The distance in bp to keep after the current bp
     * (inclusive).
     */
    public BufferedGenotypesIterator(
            VariantIterator iterator,
            int upStreamDistance,
            int downStreamDistance
    ) {

        this.iterator = iterator;
        this.upStreamDistance = upStreamDistance;
        this.downStreamDistance = downStreamDistance;

        init();

    }

    /**
     * Returns the next element. Null if iteration is finished.
     *
     * @return The next element.
     */
    public GenotypesProvider next() {

        if (currentQueue.isEmpty()) {

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

            buffer.remove(contigList.pollFirst());

            bufferSemaphore.release();

        }

        return nextElement;

    }

    /**
     * Sets up the iterator.
     */
    private void init() {

        GenotypesProvider genotypesProvider = iterator.next();

        if (genotypesProvider != null) {

            add(genotypesProvider);

        }
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

            // Load until there are enough variants in the buffer
            if (bp + upStreamDistance > currentMaxBp.get(contig)) {

                bufferSemaphore.acquire();

                int maxBp = currentMaxBp.get(contig);

                if (bp + upStreamDistance > maxBp) {

                    while (bp + LOADING_FACTOR * upStreamDistance >= maxBp) {

                        GenotypesProvider genotypesProvider = iterator.next();

                        if (genotypesProvider == null) {

                            return;

                        }

                        add(genotypesProvider);

                    }
                }

                bufferSemaphore.release();

            }

            // Remove variants outside range
            if (bp - LOADING_FACTOR * downStreamDistance > currentMaxBp.get(contig)) {

                bufferSemaphore.acquire();
                int minBp = currentMaxBp.get(contig);

                if (bp - LOADING_FACTOR * downStreamDistance > minBp) {

                    ConcurrentHashMap<Integer, ArrayList<GenotypesProvider>> contigMap = buffer.get(contig);

                    contigMap.keySet().stream()
                            .filter(
                                    tempBp -> tempBp - downStreamDistance > minBp
                            )
                            .forEach(
                                    tempBp -> contigMap.remove(tempBp)
                            );

                    currentMinBp.put(contig, minBp);

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
