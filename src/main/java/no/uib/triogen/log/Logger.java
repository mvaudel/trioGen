package no.uib.triogen.log;

import java.io.File;
import java.time.Instant;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.flat.SimpleFileWriter;

/**
 * Simple logger.
 *
 * @author Marc Vaudel
 */
public class Logger {

    private final SimpleFileWriter variantWriter;
    private final SimpleFileWriter logWriter;

    public Logger(
            File logFile
    ) {
        this(logFile, null);
    }

    public Logger(
            File logFile,
            File variantFile
    ) {

        logWriter = new SimpleFileWriter(logFile, true);
        logWriter.writeLine(
                "time",
                "type",
                "log"
        );

        if (variantFile != null) {

            variantWriter = new SimpleFileWriter(variantFile, true);
            logWriter.writeLine(
                    "time",
                    "variant",
                    "log"
            );

        } else {

            variantWriter = null;

        }
    }

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

    public void logError(String message) {

        String now = Instant.now().toString();

        logWriter.writeLine(
                now,
                "Error",
                "\"" + message.replace(IoUtils.lineSeparator, " ") + "\""
        );

    }

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
}
