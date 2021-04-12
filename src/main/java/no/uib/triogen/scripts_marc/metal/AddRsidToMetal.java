package no.uib.triogen.scripts_marc.metal;

import java.io.File;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
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

    private static final String missing = ".";

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String vcfFilePath = "/mnt/archive/utils/dbsnp/human_9606_b151_GRCh37p13/All_20180423.vcf.gz";
        String metaFilePathPattern = "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1.tbl.gz";
        String resultFilePathPattern = "/mnt/work/marc/moba/pwbw/prs/meta/{geno}/{geno}_prs1_rsid.tbl.gz";

//        String[] genos = new String[]{"child", "mother", "father"};
        String[] genos = new String[]{"mother"};

        System.out.println(Instant.now() + "    Loading variants.");

        Instant begin = Instant.now();

        // Load snp ids from metal
        HashMap<String, String> variantIdMap = new HashMap<>(0);

        for (String geno : genos) {

            String metaFilePath = metaFilePathPattern.replace("{geno}", geno);

            File metaFile = new File(metaFilePath);

            try (SimpleFileReader reader = SimpleFileReader.getFileReader(metaFile)) {

                String line = reader.readLine();

                while ((line = reader.readLine()) != null) {

                    int index = line.indexOf("\t");
                    String id = line.substring(0, index);

                    variantIdMap.put(id, null);

                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Variants lodaed, " + variantIdMap.size() + " variants to map (" + durationSeconds + " s)");

        System.out.println(Instant.now() + "    Mapping variants in " + vcfFilePath + ".");

        begin = Instant.now();

        int found = 0;

        try (VcfIterator vcfIterator = new VcfIterator(new File(vcfFilePath))) {

            VcfVariant vcfVariant;
            while ((vcfVariant = vcfIterator.next()) != null) {

                VariantInformation variantInformation = vcfVariant.getVariantInformation();

                String chr = variantInformation.contig;
                int pos = variantInformation.position;
                String rsid = variantInformation.rsId;
                String[] alleles = variantInformation.alleles;

                if (alleles.length >= 2) {

                    String[] orderedAlleles = Arrays.stream(alleles)
                            .map(
                                    allele -> allele.toUpperCase()
                            )
                            .sorted()
                            .toArray(
                                    String[]::new
                            );

                    String baseId = new StringBuilder()
                            .append(chr)
                            .append(':')
                            .append(pos)
                            .toString();

                    for (int i = 0; i < orderedAlleles.length - 1; i++) {

                        for (int j = 1; j < orderedAlleles.length; j++) {

                            String metalId = new StringBuilder(baseId)
                                    .append('_')
                                    .append(orderedAlleles[i])
                                    .append('_')
                                    .append(orderedAlleles[j])
                                    .toString();

                            if (variantIdMap.containsKey(metalId)) {

                                variantIdMap.put(metalId, rsid);

                                found++;

                            }
                        }
                    }
                }
            }
        }

        end = Instant.now();

        durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Loaded variant rsids, " + found + " variants found (" + durationSeconds + " s)");

        begin = Instant.now();

        Arrays.stream(genos)
                .parallel()
                .forEach(
                        geno -> processFile(
                                metaFilePathPattern,
                                resultFilePathPattern,
                                geno,
                                variantIdMap
                        )
                );

        end = Instant.now();

        durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Mapped ids finished (" + durationSeconds + " s)");

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
                int index = line.indexOf("\t");
                String id = line.substring(0, index);
                String rest = line.substring(index);
                writer.writeLine(id, "RSID", rest);

                while ((line = reader.readLine()) != null) {

                    index = line.indexOf("\t");
                    id = line.substring(0, index);
                    rest = line.substring(index);

                    String rsid = variantIdMap.get(id);

                    if (rsid == null) {

                        rsid = missing;

                    }

                    writer.writeLine(id, rsid, rest);

                }
            }
        }
    }
}
