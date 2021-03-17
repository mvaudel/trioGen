package no.uib.triogen.scripts_marc;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import no.uib.cell_rk.utils.SimpleFileWriter;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.genotypes.vcf.iterators.VcfIterator;
import no.uib.triogen.io.genotypes.vcf.reader.VcfVariant;
import no.uib.triogen.model.genome.VariantInformation;

/**
 * Attaches an rsid to metal results.
 *
 * @author Marc Vaudel
 */
public class AddRsidToMetal {

    private static final int batchSize = 10000000;

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String vcfFilePath = "/mnt/archive/utils/dbsnp/human_9606_b151_GRCh37p13/All_20180423.vcf.gz";
        String metaFilePathPattern = "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1.tbl.gz";
        String resultFilePathPattern = "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1_rsid.tbl.gz";

        String[] genos = new String[]{"child", "mother", "father"};

        System.out.println(Instant.now() + "    Setting up files.");

        Instant begin = Instant.now();

        String fileInPathPattern = "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1_rsid.tbl_temp.gz";

        // Process metal results
        String tempPath = fileInPathPattern;
        Arrays.stream(genos)
                .parallel()
                .forEach(
                        geno -> setupFile(
                                metaFilePathPattern,
                                tempPath,
                                geno
                        )
                );

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        int batch = 0;

        System.out.println(Instant.now() + "    Files setup finished (" + durationSeconds + " s)");

        try (VcfIterator vcfIterator = new VcfIterator(new File(vcfFilePath))) {

            VcfVariant vcfVariant;
            while ((vcfVariant = vcfIterator.next()) != null) {

                batch++;

                System.out.println(Instant.now() + "    Mapping variants for batch " + batch + ".");

                begin = Instant.now();

                HashMap<String, String> variantIdMap = new HashMap<>();

                while (vcfVariant != null && variantIdMap.size() <= batchSize) {

                    VariantInformation variantInformation = vcfVariant.getVariantInformation();

                    String chr = variantInformation.contig;
                    int pos = variantInformation.position;
                    String rsid = variantInformation.rsId;
                    String[] alleles = variantInformation.alleles;

                    TreeSet<String> orderedAlleles = Arrays.stream(alleles)
                            .map(
                                    allele -> allele.toUpperCase()
                            )
                            .collect(
                                    Collectors.toCollection(
                                            TreeSet<String>::new
                                    )
                            );

                    StringBuilder sb = new StringBuilder()
                            .append(chr)
                            .append(':')
                            .append(pos);

                    for (String allele : orderedAlleles) {

                        sb
                                .append('_')
                                .append(allele);

                    }

                    String metalId = sb.toString();

                    variantIdMap.put(metalId, rsid);

                    vcfVariant = vcfIterator.next();

                }

                end = Instant.now();

                durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                System.out.println(Instant.now() + "    Loaded variant coordinates for for batch " + batch + ", " + variantIdMap.size() + " variants (" + durationSeconds + " s)");

                System.out.println(Instant.now() + "    Mapping ids for batch " + batch + ".");

                begin = Instant.now();

                String inPath = fileInPathPattern;
                String fileOutPathPattern = vcfVariant == null ? resultFilePathPattern : "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1_rsid.tbl_temp" + batch + ".gz";

                Arrays.stream(genos)
                        .parallel()
                        .forEach(
                                geno -> processFile(
                                        inPath,
                                        fileOutPathPattern,
                                        geno,
                                        variantIdMap
                                )
                        );

                end = Instant.now();

                durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                System.out.println(Instant.now() + "    Mapped ids for for batch " + batch + ", " + variantIdMap.size() + " variants (" + durationSeconds + " s)");

                System.out.println(Instant.now() + "    Deleting temp file for batch " + batch + ".");

                begin = Instant.now();

                Arrays.stream(genos)
                        .parallel()
                        .map(
                                geno -> new File(inPath.replace("{geno}", geno))
                        )
                        .forEach(
                                file -> file.delete()
                        );

                fileInPathPattern = fileOutPathPattern;

                end = Instant.now();

                durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

                System.out.println(Instant.now() + "    Temp file deleted for batch " + batch + " (" + durationSeconds + " s)");

            }
        }
    }

    /**
     * Processes a file.
     *
     * @param tempFilePathPattern The temp file name pattern.
     * @param resultFilePathPattern The result file name pattern.
     * @param geno The geno key.
     * @param variantIdMap The variant id to rsid map.
     */
    private static void processFile(
            String tempFilePathPattern,
            String resultFilePathPattern,
            String geno,
            HashMap<String, String> variantIdMap
    ) {

        String tempFilePath = tempFilePathPattern.replace("{geno}", geno);
        String resultFilePath = resultFilePathPattern.replace("{geno}", geno);

        File tempFile = new File(tempFilePath);
        File resultFile = new File(resultFilePath);

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(tempFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(resultFile, true)) {

                String line = reader.readLine();
                writer.writeLine(line);

                while ((line = reader.readLine()) != null) {

                    int index = line.indexOf("\t");
                    String id = line.substring(0, index);

                    String rsid = variantIdMap.get(id);

                    if (rsid != null) {

                        String rest = line.substring(index);
                        index = rest.indexOf("\t");
                        rest = line.substring(index);
                        writer.writeLine(id, rsid, rest);

                    } else {

                        writer.writeLine(line);

                    }
                }
            }
        }
    }

    /**
     * Processes a file.
     *
     * @param metaFilePathPattern The meta file name pattern.
     * @param resultFilePathPattern The result file name pattern.
     * @param geno The geno key.
     * @param variantIdMap The variant id to rsid map.
     */
    private static void setupFile(
            String metaFilePathPattern,
            String resultFilePathPattern,
            String geno
    ) {

        String metaFilePath = metaFilePathPattern.replace("{geno}", geno);
        String resultFilePath = resultFilePathPattern.replace("{geno}", geno);

        File metaFile = new File(metaFilePath);
        File resultFile = new File(resultFilePath);

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(metaFile)) {

            try (SimpleFileWriter writer = new SimpleFileWriter(resultFile, true)) {

                String line = reader.readLine();

                int index = line.indexOf("\t");
                String id = line.substring(0, index);
                String rest = line.substring(index);
                writer.writeLine(id, "RSID", rest);

                while ((line = reader.readLine()) != null) {

                    index = line.indexOf("\t");
                    id = line.substring(0, index);
                    rest = line.substring(index);

                    String rsid = ".";

                    writer.writeLine(id, rsid, rest);

                }
            }
        }
    }
}
