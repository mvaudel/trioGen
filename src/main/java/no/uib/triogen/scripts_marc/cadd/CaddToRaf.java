package no.uib.triogen.scripts_marc.cadd;

import java.io.File;
import java.io.RandomAccessFile;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.time.Instant;
import java.util.ArrayList;
import no.uib.triogen.io.IoUtils;
import static no.uib.triogen.io.IoUtils.ENCODING;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.utils.CompressionUtils;
import no.uib.triogen.utils.TempByteArray;

/**
 * Puts the cadd database in a raf file.
 *
 * java -Xmx32G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.scripts_marc.cadd.CaddToRaf
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

        File caddFile = new File("/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz");

        String outputPath = "/mnt/work2/utils/cadd/GRCh37/whole_genome_SNVs_inclAnno";

        try {

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(caddFile, false)) {

                String line = reader.readLine();
                String header = reader.readLine();

                String lastChromosome = "-1";
                int lastBp = 0;

                ArrayList<byte[]> keys = new ArrayList<>();
                ArrayList<Long> indexes = new ArrayList<>();
                int keyMapSize = 0;

                while ((line = reader.readLine()) != null) {

                    String[] lineSplit = line.split("\t");

                    String chromosome = lineSplit[0];
                    String bpString = lineSplit[1];
                    String ref = lineSplit[2];
                    String alt = lineSplit[3];

                    int bp = Integer.parseInt(bpString);

                    if (bp - lastBp > 100000) {

                        System.out.println(Instant.now() + " -    " + bp);

                        lastBp = bp;

                    }

                    String key = String.join("_", chromosome, bpString, ref, alt);

                    if (!lastChromosome.equals(chromosome)) {

                        if (raf != null) {

                            ByteBuffer byteBuffer = ByteBuffer.allocate(keyMapSize);

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
                    keyMapSize += Integer.BYTES + keyBytes.length + Long.BYTES;

                    byte[] lineBytes = line.getBytes(IoUtils.ENCODING);
                    TempByteArray array = CompressionUtils.zstdCompress(lineBytes);

                    raf.write(lineBytes.length);
                    raf.write(array.array, 0, array.length);

                }
            }

        } catch (Throwable t) {

            t.printStackTrace();

        }
    }

}
