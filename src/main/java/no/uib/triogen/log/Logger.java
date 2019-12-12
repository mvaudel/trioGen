package no.uib.triogen.log;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.TrioGen;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Simple logger.
 *
 * @author Marc Vaudel
 */
public class Logger {

    /**
     * Writer for the variant log.
     */
    private final SimpleFileWriter variantWriter;
    /**
     * Writer for the general log.
     */
    private final SimpleFileWriter logWriter;

    /**
     * Constructor.
     *
     * @param logFile the file where to write the log
     */
    public Logger(
            File logFile
    ) {
        this(logFile, null);
    }

    /**
     * Constructor.
     *
     * @param logFile the file where to write the log
     * @param variantFile the file where to write the variants log
     */
    public Logger(
            File logFile,
            File variantFile
    ) {

        logWriter = new SimpleFileWriter(logFile, true);
        logWriter.writeLine("# TrioGen version: " + TrioGen.getVersion());
        logWriter.writeLine(
                "time",
                "type",
                "log"
        );

        if (variantFile != null) {

            variantWriter = new SimpleFileWriter(variantFile, true);
            variantWriter.writeLine("# TrioGen version " + TrioGen.getVersion());
            variantWriter.writeLine(
                    "time",
                    "variant",
                    "log"
            );

        } else {

            variantWriter = null;

        }
    }

    /**
     * Logs a message.
     *
     * @param message the message
     */
    public void logMessage(
            String message
    ) {

        String now = Instant.now().toString();

        logWriter.writeLine(
                now,
                "Message",
                "\"" + message.replace(IoUtils.lineSeparator, " ") + "\""
        );

        System.out.println(now + " - " + message);

    }

    /**
     * Logs an error.
     *
     * @param message the error message
     */
    public void logError(String message) {

        String now = Instant.now().toString();

        logWriter.writeLine(
                now,
                "Error",
                "\"" + message.replace(IoUtils.lineSeparator, " ") + "\""
        );

    }

    /**
     * Log for a variant.
     *
     * @param variantId the id of the variant
     * @param message the message
     */
    public void logVariant(
            String variantId,
            String message
    ) {

        if (variantWriter != null) {

            String now = Instant.now().toString();

            variantWriter.writeLine(
                    now,
                    variantId,
                    message
            );

        }
    }

    /**
     * Closes the logger.
     */
    public void close() {

        logWriter.close();

        if (variantWriter != null) {

            variantWriter.close();

        }
    }
}
