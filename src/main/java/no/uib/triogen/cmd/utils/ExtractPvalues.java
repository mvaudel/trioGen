package no.uib.triogen.cmd.utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileReader;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Extracts the p-values by pheno from a transmission result file. For personal
 * use only.
 *
 * @author Marc Vaudel
 */
public class ExtractPvalues {

    /**
     * Main method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        (new ExtractPvalues()).run(args[0], Arrays.copyOfRange(args, 1, args.length));

    }

    private void run(String resultsFilePath, String[] phenos) {

        String stem = resultsFilePath.endsWith(".gz") ? resultsFilePath.substring(0, resultsFilePath.indexOf(".gz")) : resultsFilePath;

        File resultFile = new File(resultsFilePath);

        HashMap<String, SimpleFileWriter> fileWriters = Arrays.stream(phenos)
                .collect(
                        Collectors.toMap(
                                pheno -> pheno,
                                pheno -> new SimpleFileWriter(
                                        new File(stem + "_" + pheno + ".gz"),
                                        true
                                ),
                                (a, b) -> a,
                                HashMap::new)
                );
        
        try (SimpleFileReader reader = SimpleFileReader.getFileReader(resultFile)) {
            
            String line = reader.readLine();
            String[] headerSplit = line.split(IoUtils.separator);
            
            ArrayList<Integer> columns = new ArrayList<>();
            StringBuilder stringBuilder = new StringBuilder();
            
            for (int i = 0 ; i < headerSplit.length ; i++) {
                
                String column = headerSplit[i];
                
                if (column.endsWith("_p") || i < 5) {
                    
                    columns.add(i);
                    
                }
                
                if (i > 0) {
                    
                stringBuilder.append(IoUtils.separator);
                    
                }
                
                String newColum = column.replaceAll("β", "B");
                
                stringBuilder.append(newColum);
            }
            
            String newHeader = stringBuilder.toString();
            
            for (SimpleFileWriter writer : fileWriters.values()) {
                
                writer.writeLine(newHeader);
                
            }
            
            while((line = reader.readLine()) != null) {
                
                String[] lineSplit = line.split(IoUtils.separator);
                
                String pheno = lineSplit[0];
                
                SimpleFileWriter writer = fileWriters.get(pheno);
                
                if (writer != null) {
                    
                    String newLine = columns.stream()
                            .map(
                                    i -> lineSplit[i]
                            )
                            .collect(
                                    Collectors.joining(IoUtils.separator)
                            );
                    
                    writer.writeLine(newLine);
                    
                }                
            }            
        }
        
        fileWriters.values()
                .forEach(
                        writer -> writer.close()
                );

    }

}
