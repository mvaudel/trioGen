package no.uib.triogen.io.flat.mapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.ld.LdMatrixWriter;
import no.uib.triogen.utils.SimpleSemaphore;

/**
 * This class maps a file to memory using evenly distributed blocks of the given
 * block size. Threads accessing the file block by block.
 *
 * @author Marc Vaudel
 */
public class MemoryMappedFile implements AutoCloseable {

    /**
     * The size to use for the blocks.
     */
    private final long blockSize;
    /**
     * The random access file.
     */
    private final RandomAccessFile raf;
    /**
     * The channel to the file.
     */
    private final FileChannel fc;
    /**
     * Semaphore to synchronize threads.
     */
    private final SimpleSemaphore[] bufferSemaphores;
    /**
     * The mapped byte buffer.
     */
    private final MappedByteBuffer[] mappedByteBuffers;

    /**
     * Constructor.
     *
     * @param file The file to map.
     * @param offset The start index of the section of the file to map.
     * @param length The length of the section of the file to map.
     * @param blockSize The block size to use in number of bytes.
     *
     * @throws IOException Exception thrown if an error occurred while reading
     * the file.
     */
    public MemoryMappedFile(
            File file,
            long offset,
            long length,
            long blockSize
    ) throws IOException {

        this.blockSize = blockSize;

        raf = new RandomAccessFile(file, "r");

        try {

            fc = raf.getChannel();

            int nBuffers = (int) (length / blockSize);
            int rest = (int) (length - nBuffers * blockSize);

            if (rest > 0) {

                nBuffers++;

            }

            bufferSemaphores = new SimpleSemaphore[nBuffers];
            mappedByteBuffers = new MappedByteBuffer[nBuffers];

            for (int i = 0; i < nBuffers; i++) {

                bufferSemaphores[i] = new SimpleSemaphore(1);

                long begin = i * blockSize + LdMatrixWriter.HEADER_LENGTH;
                long bufferSizeI = i < nBuffers - 1 || rest == 0 ? blockSize : rest;

                mappedByteBuffers[i] = fc.map(
                        FileChannel.MapMode.READ_ONLY,
                        begin,
                        bufferSizeI
                );
            }

        } finally {

            raf.close();

        }
    }

    /**
     * Returns a buffer at the given position in the section of the file mapped.
     * 0 is the beginning of the mapped section.
     *
     * @param index The index at which to start the buffer.
     *
     * @return A buffer at the given position in the section of the file mapped.
     */
    public MiniBuffer getBuffer(long index) {

        return new MiniBuffer(index);

    }

    @Override
    public void close() throws Exception {

        fc.close();

        for (MappedByteBuffer mappedByteBuffer : mappedByteBuffers) {

            IoUtils.closeBuffer(mappedByteBuffer);

        }
    }

    /**
     * Buffer iterating over the blocks of the file.
     */
    public class MiniBuffer implements AutoCloseable {

        /**
         * The index of the block being buffered.
         */
        private int blockIndex;
        /**
         * Mapped byte buffer being buffered.
         */
        private MappedByteBuffer mappedByteBuffer;
        /**
         * Placeholder for the array of bytes used to decode overlapping
         * primitive data types.
         */
        private final byte[] tempArray = new byte[Double.BYTES];
        /**
         * Placeholder for the byte buffer used to decode overlapping primitive
         * data types.
         */
        private final ByteBuffer byteBuffer = ByteBuffer.allocate(Double.BYTES);

        /**
         * Constructor.
         *
         * @param index The index in the section of the file mapped to buffer.
         */
        private MiniBuffer(
                long index
        ) {

            blockIndex = (int) (index / blockSize);
            int indexInBuffer = (int) (index - blockIndex * blockSize);

            mappedByteBuffer = mappedByteBuffers[blockIndex];
            SimpleSemaphore bufferSemaphore = bufferSemaphores[blockIndex];

            bufferSemaphore.acquire();

            ((Buffer) mappedByteBuffer).position(indexInBuffer);

        }

        /**
         * Returns the integer at the current position in the file.
         *
         * @return The integer at the current position in the file.
         */
        public int getInt() {

            int remaining = mappedByteBuffer.remaining();

            if (remaining == 0) {

                nextBlock();
                remaining = mappedByteBuffer.remaining();

            }

            if (remaining >= Integer.BYTES) {

                return mappedByteBuffer.getInt();

            } else {

                getOverlappingBytes(tempArray, Integer.BYTES);

                byteBuffer.put(tempArray, 0, Integer.BYTES);

                int result = byteBuffer.getInt();

                byteBuffer.position(0);

                return result;

            }
        }

        /**
         * Returns the long at the current position in the file.
         *
         * @return The long at the current position in the file.
         */
        public long getLong() {

            int remaining = mappedByteBuffer.remaining();

            if (remaining == 0) {

                nextBlock();
                remaining = mappedByteBuffer.remaining();

            }

            if (remaining >= Long.BYTES) {

                return mappedByteBuffer.getLong();

            } else {

                getOverlappingBytes(tempArray, Long.BYTES);

                byteBuffer.put(tempArray, 0, Long.BYTES);

                long result = byteBuffer.getLong();

                byteBuffer.position(0);

                return result;

            }
        }

        /**
         * Returns the float at the current position in the file.
         *
         * @return The float at the current position in the file.
         */
        public float getFloat() {

            int remaining = mappedByteBuffer.remaining();

            if (remaining == 0) {

                nextBlock();
                remaining = mappedByteBuffer.remaining();

            }

            if (remaining >= Float.BYTES) {

                return mappedByteBuffer.getFloat();

            } else {

                getOverlappingBytes(tempArray, Float.BYTES);

                byteBuffer.put(tempArray, 0, Float.BYTES);

                float result = byteBuffer.getFloat();

                byteBuffer.position(0);

                return result;

            }
        }

        /**
         * Returns the double at the current position in the file.
         *
         * @return The double at the current position in the file.
         */
        public double getDouble() {

            int remaining = mappedByteBuffer.remaining();

            if (remaining == 0) {

                nextBlock();
                remaining = mappedByteBuffer.remaining();

            }

            if (remaining >= Double.BYTES) {

                return mappedByteBuffer.getDouble();

            } else {

                getOverlappingBytes(tempArray, Double.BYTES);

                byteBuffer.put(tempArray, 0, Double.BYTES);

                double result = byteBuffer.getDouble();

                byteBuffer.position(0);

                return result;

            }
        }

        /**
         * Maps the bytes at the current position to the given array.
         *
         * @param bytes The byte array to fill.
         *
         * @return The current buffer.
         */
        public MiniBuffer get(byte[] bytes) {

            int remaining = mappedByteBuffer.remaining();

            if (remaining == 0) {

                nextBlock();
                remaining = mappedByteBuffer.remaining();

            }

            if (remaining >= bytes.length) {

                mappedByteBuffer.get(bytes);

            } else {

                getOverlappingBytes(bytes, bytes.length);

            }

            return this;

        }

        private void getOverlappingBytes(
                byte[] buffer,
                int nBytes
        ) {

            int loaded = mappedByteBuffer.remaining();

            mappedByteBuffer.get(buffer, 0, loaded);

            while (loaded < nBytes) {

                nextBlock();

                int length = Math.min(nBytes - loaded, mappedByteBuffer.remaining());

                mappedByteBuffer.get(buffer, loaded, length);

                loaded += length;

            }
        }

        /**
         * Moves to the next block.
         */
        private void nextBlock() {

            SimpleSemaphore bufferSemaphore = bufferSemaphores[blockIndex];
            bufferSemaphore.release();

            blockIndex++;

            mappedByteBuffer = mappedByteBuffers[blockIndex];
            bufferSemaphore = bufferSemaphores[blockIndex];

            bufferSemaphore.acquire();
            ((Buffer) mappedByteBuffer).position(0);

        }

        @Override
        public void close() {

            SimpleSemaphore bufferSemaphore = bufferSemaphores[blockIndex];
            bufferSemaphore.release();

        }

    }

}
