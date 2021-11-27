/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uib.triogen.utils.cli;

/**
 * Interface for a CLI option.
 *
 * @author Marc Vaudel
 */
public interface CliOption {
    
    /**
     * Returns the short option name.
     * 
     * @return The short option name.
     */
    public String getOption();
    
    /**
     * Returns the long option name.
     * 
     * @return The long option name.
     */
    public String getLongOption();
    
}
