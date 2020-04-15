package no.uib.triogen.utils;

/**
 * Semaphore allowing delays.
 *
 * @author Marc Vaudel
 */
public class SimpleTimeSemaphore {

    /**
     * The time to delay in ms.
     */
    private int delay = 10;

    /**
     * Delays of the
     *
     * @throws InterruptedException
     */
    public synchronized void delay() throws InterruptedException {

        wait(delay);
    }

    /**
     * Increases the delay.
     */
    public void increaseDelay() {

        if (delay < 5000) {

            delay *= 1.1;

        }
    }

    /**
     * Sets the delay in ms.
     *
     * @param delay The delay in ms.
     */
    public void setDelay(
            int delay
    ) {
        this.delay = delay;
    }

}
