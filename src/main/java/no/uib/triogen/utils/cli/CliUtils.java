package no.uib.triogen.utils.cli;

import org.apache.commons.cli.CommandLine;

/**
 * Convenience methods for command line handling.
 *
 * @author Marc Vaudel
 */
public class CliUtils {
    
    /**
     * Returns a boolean indicating whether the command line has the given option.
     * 
     * @param aLine The command line.
     * @param cliOption The command line option.
     * 
     * @return A boolean indicating whether the command line has the given option.
     */
    public static boolean hasOption(
            CommandLine aLine,
            CliOption cliOption
    ) {
        
        return aLine.hasOption(cliOption.getOption()) || aLine.hasOption(cliOption.getLongOption());
        
    }
    
    /**
     * Returns the value at a given option in the given command line.
     * 
     * @param aLine The command line.
     * @param cliOption The command line option.
     * 
     * @return The value at a given option in the given command line.
     */
    public static String getOptionValue(
            CommandLine aLine,
            CliOption cliOption
    ) {
        
        if (aLine.hasOption(cliOption.getOption())) {
            
            return aLine.getOptionValue(cliOption.getOption());
            
        }
        if (aLine.hasOption(cliOption.getLongOption())) {
            
            return aLine.getOptionValue(cliOption.getLongOption());
            
        }
        
        return null;
        
    }
}
