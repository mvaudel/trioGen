package no.uib.triogen.scripts_marc.cadd;

import java.io.File;
import java.io.RandomAccessFile;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.time.Instant;
import java.util.ArrayList;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.utils.CompressionUtils;
import no.uib.triogen.utils.TempByteArray;

/**
 * Puts the cadd database in a raf file.
 *
 * java -Xmx64G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar
 * no.uib.triogen.scripts_marc.cadd.CaddToRaf
 *
 * @author Marc Vaudel
 */
public class CaddToRaf {

    /**
     * The magic number of to use to identify the supported files.
     */
    public static final byte[] MAGIC_NUMBER = getMagicNumber();
    /**
     * The length of the file header.
     */
    public static final int HEADER_LENGTH = MAGIC_NUMBER.length + Long.BYTES;
    /**
     * The random access file to write to.
     */
    private static RandomAccessFile raf = null;

    /**
     * Returns the magic number.
     *
     * @return The magic number.
     */
    public static byte[] getMagicNumber() {

        try {

            String magicName = "Triogen.cadd.1.1";
            return magicName.getBytes(ENCODING);

        } catch (UnsupportedEncodingException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File caddFile = new File("/mnt/archive/marc/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz");

        String outputPath = "/mnt/archive/marc/cadd/GRCh37/whole_genome_SNVs_inclAnno";

        try {

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(caddFile, false)) {

                String line = reader.readLine();
                String header = reader.readLine();

                String lastChromosome = "-1";
                int lastBp = 0;
                String currentKey = "";
                ArrayList<String> currentLines = new ArrayList<>();

                ArrayList<byte[]> keys = new ArrayList<>();
                ArrayList<Long> indexes = new ArrayList<>();
                int keyMapSize = 0;

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String chromosome = lineSplit[0];
                    String bpString = lineSplit[1];
                    String ref = lineSplit[2];
                    String alt = lineSplit[3];

                    String key = String.join("_", chromosome, bpString, ref, alt);

                    if (currentKey.equals("")) {

                        currentKey = key;

                    }

                    if (!key.equals(currentKey)) {

                        int bp = Integer.parseInt(bpString);

                        if (bp - lastBp > 100000) {

                            System.out.println(Instant.now() + " -    " + bp);

                            lastBp = bp;

                        }

                        if (!lastChromosome.equals(chromosome)) {

                            if (raf != null) {

                                System.out.println(Instant.now() + " - Chromosome " + lastChromosome + " done, saving index");

                                ByteBuffer byteBuffer = ByteBuffer.allocate(keyMapSize);

                                System.out.println(Instant.now() + " - keyMapSize: " + keyMapSize + " indexes: " + indexes.size());

                                for (int i = 0; i < indexes.size(); i++) {

                                    byte[] keyBytes = keys.get(i);
                                    long index = indexes.get(i);

                                    byteBuffer.putInt(keyBytes.length);
                                    byteBuffer.put(keyBytes);
                                    byteBuffer.putLong(index);

                                }

                                long footerPosition = raf.getFilePointer();

                                byte[] indexesBytes = byteBuffer.array();
                                TempByteArray array = CompressionUtils.zstdCompress(indexesBytes);

                                raf.write(indexesBytes.length);
                                raf.write(array.array, 0, array.length);

                                byte[] headerBytes = header.getBytes(IoUtils.ENCODING);
                                array = CompressionUtils.zstdCompress(headerBytes);

                                raf.write(headerBytes.length);
                                raf.write(array.array, 0, array.length);

                                raf.seek(0);
                                raf.write(MAGIC_NUMBER);
                                raf.writeLong(footerPosition);

                                return;

                            }

                            System.out.println(Instant.now() + " - Processing Chromosome " + chromosome);

                            File outputFile = new File(outputPath + "_" + chromosome);

                            raf = new RandomAccessFile(outputFile, "rw");
                            raf.seek(HEADER_LENGTH);

                            lastChromosome = chromosome;
                            lastBp = bp;

                            keys.clear();
                            indexes.clear();
                            keyMapSize = 0;

                        }

                        byte[] keyBytes = key.getBytes(IoUtils.ENCODING);
                        keys.add(keyBytes);
                        indexes.add(raf.getFilePointer());

                        int toAdd = Integer.BYTES + keyBytes.length + Long.BYTES;

                        if (((long) keyMapSize) + toAdd > Integer.MAX_VALUE) {

                            throw new IllegalArgumentException("Key buffer higher than int");

                        }

                        keyMapSize += toAdd;

                        String lines = currentLines.stream()
                                .collect(
                                        Collectors.joining(System.lineSeparator())
                                );
                        byte[] linesBytes = lines.getBytes(IoUtils.ENCODING);
                        TempByteArray array = CompressionUtils.zstdCompress(linesBytes);

                        raf.write(linesBytes.length);
                        raf.write(array.array, 0, array.length);

                        currentLines.clear();
                        currentKey = key;

                    }

                    currentLines.add(line);

                }

                byte[] keyBytes = currentKey.getBytes(IoUtils.ENCODING);
                keys.add(keyBytes);
                indexes.add(raf.getFilePointer());
                keyMapSize += Integer.BYTES + keyBytes.length + Long.BYTES;

                String lines = currentLines.stream()
                        .collect(
                                Collectors.joining(System.lineSeparator())
                        );
                byte[] linesBytes = lines.getBytes(IoUtils.ENCODING);
                TempByteArray array = CompressionUtils.zstdCompress(linesBytes);

                raf.write(linesBytes.length);
                raf.write(array.array, 0, array.length);

                System.out.println(Instant.now() + " - Chromosome " + lastChromosome + " done, saving index");

                ByteBuffer byteBuffer = ByteBuffer.allocate(keyMapSize);

                for (int i = 0; i < indexes.size(); i++) {

                    byte[] tempKeyBytes = keys.get(i);
                    long index = indexes.get(i);

                    byteBuffer.putInt(tempKeyBytes.length);
                    byteBuffer.put(tempKeyBytes);
                    byteBuffer.putLong(index);

                }

                long footerPosition = raf.getFilePointer();

                byte[] indexesBytes = byteBuffer.array();
                array = CompressionUtils.zstdCompress(indexesBytes);

                raf.write(indexesBytes.length);
                raf.write(array.array, 0, array.length);

                byte[] headerBytes = header.getBytes(IoUtils.ENCODING);
                array = CompressionUtils.zstdCompress(headerBytes);

                raf.write(headerBytes.length);
                raf.write(array.array, 0, array.length);

                raf.seek(0);
                raf.write(MAGIC_NUMBER);
                raf.writeLong(footerPosition);

            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }

}
