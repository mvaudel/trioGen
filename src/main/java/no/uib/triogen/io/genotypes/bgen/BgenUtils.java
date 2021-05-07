package no.uib.triogen.io.genotypes.bgen;

import java.io.UnsupportedEncodingException;
import java.util.concurrent.ConcurrentHashMap;
import static no.uib.triogen.io.IoUtils.ENCODING;

/**
 * Utils for the reading of bgen files.
 *
 * @author Marc Vaudel
 */
public class BgenUtils {
    
    /**
     * The magic number of to use to identify the supported files.
     */
    public static final byte[] MAGIC_NUMBER = getMagicNumber();
    /**
     * The identifier for the TrioGen formatted files.
     */
    public static final byte[] IDENTIFIER = getIdentifier();
    /**
     * Cache for powers of two used in probability encoding.
     */
    private static final ConcurrentHashMap<Integer, int[]> powersCache = new ConcurrentHashMap<>(1);
    
    /**
     * Returns the powers of two for probability encoding in an array.
     * 
     * @param nBits The number of bits used for encoding.
     * 
     * @return The powers of two for probability encoding in an array.
     */
    public static int[] getPowers(
    int nBits
    ) {
        
        int[] powers = powersCache.get(nBits);
        
        if (powers == null) {
            
            powers = new int[nBits];

                int value = 1;

                for (int i = 0; i < nBits; i++) {

                    powers[i] = value;

                    value *= 2;

                }
                
                powersCache.put(nBits, powers);
            
        }
        
        return powers;
        
    }

    /**
     * Returns the magic number.
     *
     * @return The magic number.
     */
    public static byte[] getMagicNumber() {

        try {

            String magicName = "bgen";
            byte[] magicNumber = magicName.getBytes(ENCODING);
            
            if (magicNumber.length != 4) {
                
                throw new IllegalArgumentException("Bgen magic number should be of length 4.");
                
            }
            
            return magicNumber;

        } catch (UnsupportedEncodingException e) {

            throw new RuntimeException(e);

        }
    }

    /**
     * Returns the identifier for TrioGen-generated files.
     *
     * @return The identifier for TrioGen-generated files.
     */
    public static byte[] getIdentifier() {

        try {

            String idString = "TrioGen-0.5.0";
            byte[] idBytes = idString.getBytes(ENCODING);
            
            return idBytes;

        } catch (UnsupportedEncodingException e) {

            throw new RuntimeException(e);

        }
    }
    
    public static boolean checkMagicNumber(byte[] fileMagicNumber) {
        
        if (fileMagicNumber.length != 4) {
            
            return false;
            
        }
        
        // Backward compatibility
        
        boolean allZero = true;
        
        for (int i = 0 ; i < 3 ; i++) {
            
            if (fileMagicNumber[i] != 0b0) {
                
                allZero = false;
                break;
                
            }
        }
        
        if (allZero) {
            
            return true;
            
        }
        
        // Check bgen
        
        for (int i = 0 ; i < 3 ; i++) {
            
            if (fileMagicNumber[i] != MAGIC_NUMBER[i]) {
                
                return false;
                
            }
        }
        
        return true;
        
    }
}
