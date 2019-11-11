package no.uib.triogen.model.pheno;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import no.uib.triogen.io.Utils;
import no.uib.triogen.io.flat.SimpleFileReader;

/**
 * This class parses phenotypes from a pheno file.
 *
 * @author Marc Vaudel
 */
public class PhenotypesHandler {

    /**
     * Name of the column containing the child ids.
     */
    public final static String childIdColumn = "childId";

    /**
     * Phenotype name to child id to phenotype value map.
     */
    public final HashMap<String, HashMap<String, Double>> phenoMap = new HashMap<>();
    
    /**
     * The number of children for which a phenotype was found.
     */
    public final int nChildren;

    /**
     * Constructor from a phenotypes file.
     * 
     * @param phenoFile a file containing the phenotypes
     * @param phenoNames the names of the phenotypes to parse
     */
    public PhenotypesHandler(
            File phenoFile,
            String[] phenoNames
    ) {
        
        for (String phenoName : phenoNames) {
            
            phenoMap.put(phenoName, new HashMap<>(1000));
            
        }

        try (SimpleFileReader phenoReader = SimpleFileReader.getFileReader(phenoFile)) {

            String line = phenoReader.readLine();

            if (line == null) {

                throw new IllegalArgumentException(
                        "Pheno file is empty."
                );

            }

            String[] lineContent = line.split(Utils.separator);
            
            int nColumns = lineContent.length;

            if (nColumns < 2) {

                throw new IllegalArgumentException(
                        "Only one column found in pheno file. Please verify that the pheno file is correctly formatted."
                );

            }

            HashSet<String> phenoNamesSet = Arrays.stream(phenoNames)
                    .collect(Collectors.toCollection(HashSet::new));
            
            int idColumnIndex = -1;
            HashMap<String, Integer> phenoColumnIndexMap = new HashMap<>(phenoNames.length);

            for (int i = 0; i < lineContent.length; i++) {
                
                String cellContent = lineContent[i];

                if (cellContent.equals(PhenotypesHandler.childIdColumn)) {

                    idColumnIndex = i;

                } else if (phenoNamesSet.contains(cellContent)) {
                    
                    phenoColumnIndexMap.put(cellContent, i);
                    
                }
            }

            if (idColumnIndex == -1) {

                throw new IllegalArgumentException(
                        "The phenotypes file does not contain a \"" + PhenotypesHandler.childIdColumn + "\" column."
                );

            }

            String[] missingPhenos = Arrays.stream(phenoNames)
                    .filter(
                            phenoName -> !phenoColumnIndexMap.containsKey(phenoName)
                    )
                    .toArray(String[]::new);

            if (missingPhenos.length > 0) {

                throw new IllegalArgumentException(
                        "The following phenotypes were not found in the phenotypes file: " + String.join(", ", missingPhenos) + "."
                );

            }
            
            int lineNumber = 1;
            
            while ((line = phenoReader.readLine()) != null) {
                
                lineNumber++;
                
                lineContent = line.split(Utils.separator);
                
                if (lineContent.length != nColumns) {
                    
                    throw new IllegalArgumentException(
                        "Unexpected number of columns at line " + lineNumber + " (Expected: " + nColumns + ", Found: " + lineContent.length + ")."
                );
                    
                }
                
                String childId = lineContent[idColumnIndex];
                
                for (Entry<String, Integer> phenoColumn : phenoColumnIndexMap.entrySet()) {
                    
                    String valueString = lineContent[phenoColumn.getValue()];
                    
                    try {
                        
                        double value = Double.parseDouble(valueString);
                        
                        phenoMap.get(phenoColumn.getKey()).put(childId, value);
                        
                    } catch (Exception e) {
                        
                        throw new IllegalArgumentException(
                        "The value for phenotype " + phenoColumn.getKey() + " at line " + lineNumber + " (" + valueString + ") could not be parsed as a number."
                );
                        
                    }
                }
            }
            
            nChildren = lineNumber - 1;
            
        }
    }
}
