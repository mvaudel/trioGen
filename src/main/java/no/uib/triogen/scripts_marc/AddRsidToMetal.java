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

        // Load variant coordinates
        System.out.println(Instant.now() + "    Loading variant coordinates from '" + vcfFilePath + "' .");

        Instant begin = Instant.now();

        HashMap<String, String> variantIdMap = new HashMap<>();

        try (VcfIterator vcfIterator = new VcfIterator(new File(vcfFilePath))) {

            VcfVariant vcfVariant;
            while ((vcfVariant = vcfIterator.next()) != null) {

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

            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    Loaded variant coordinates for " + variantIdMap.size() + " variants (" + durationSeconds + " s)");

        // Process metal results
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
    }

    /**
     * Processes a file.
     * 
     * @param metaFilePathPattern The meta file name pattern.
     * @param resultFilePathPattern The result file name pattern.
     * @param geno The geno key.
     * @param variantIdMap The variant id to rsid map.
     */
    private static void processFile(
            String metaFilePathPattern,
            String resultFilePathPattern,
            String geno,
            HashMap<String, String> variantIdMap
    ) {

        Instant begin = Instant.now();

        String metaFilePath = metaFilePathPattern.replace("{geno}", geno);
        String resultFilePath = resultFilePathPattern.replace("{geno}", geno);

        System.out.println(Instant.now() + "    Processing '" + metaFilePath + "' .");

        File metaFile = new File(metaFilePath);
        File resultFile = new File(resultFilePath);
        
        int found = 0;

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

                    String rsid = variantIdMap.get(id);

                    if (rsid == null) {

                        rsid = ".";

                    } else {
                        
                        found++;
                        
                    }

                    writer.writeLine(id, rsid, rest);

                }
            }
        }

        Instant end = Instant.now();

        long durationSeconds = end.getEpochSecond() - begin.getEpochSecond();

        System.out.println(Instant.now() + "    " + found + " rsids mapped in  '" + metaFilePath + "' (" + durationSeconds + " s)");

    }
}
