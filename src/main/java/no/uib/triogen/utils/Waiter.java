package no.uib.triogen.utils;

/**
 * A class to wait.
 *
 * @author Marc Vaudel
 */
public class Waiter {
    
    public synchronized void delay(long time) throws InterruptedException {
        
        this.wait(time);
        
    }

}
