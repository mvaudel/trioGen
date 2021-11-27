package no.uib.triogen.cmd.prs_prune;

import java.io.File;
import no.uib.triogen.cmd.ld_pruning.LdPruningOptions;
import no.uib.triogen.cmd.prs_score.PrsScoreOptions;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsPruner;
import no.uib.triogen.utils.cli.CliUtils;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class PrsPruneOptionsBean {

    /**
     * The training file.
     */
    public final File trainingFile;
    /**
     * The ld matrix file path.
     */
    public final String ldMatrixFilePath;
    /**
     * Name of the variant identifier column in the training file.
     */
    public String snpIdColumn = "variantId";
    /**
     * Name of the contig column in the training file.
     */
    public String chrColumn = "contig";
    /**
     * Name of the position column in the training file.
     */
    public String posColumn = "position";
    /**
     * Name of the reference allele column in the training file.
     */
    public String refColumn = "otherAllele";
    /**
     * Name of the effect allele column in the training file.
     */
    public String eaColumn = "testedAllele";
    /**
     * Names of the variables.
     */
    public String[] variables = new String[]{"c", "m", "f"};
    /**
     * Pattern for the effect size column.
     */
    public String betaPattern = Model.cmf.name() + ".B" + PrsPruner.VARIABLE_WILDCARD;
    /**
     * Pattern for the standard error column.
     */
    public String sePattern = Model.cmf.name() + ".B.se" + PrsPruner.VARIABLE_WILDCARD;
    /**
     * Pattern for the p-value column.
     */
    public String pPattern = Model.cmf.name() + ".B.p" + PrsPruner.VARIABLE_WILDCARD;
    /**
     * The file where to write the output.
     */
    public final File destinationFile;
    /**
     * LD R2 value over which two hits cannot be considered independent.
     */
    public double ldLocusThreshold = 0.05;
    /**
     * LD R2 value over which two hits are considered identical.
     */
    public double ldTopHitThreshold = 0.8;
    /**
     * The minimal number of variants required for a locus.
     */
    public int nSnpPerLocusThreshold = 5;
    /**
     * The highest p-value to consider.
     */
    public double pValueThreshold = 0.05;
    /**
     * The allele frequency threshold.
     */
    public double afThreshold = Double.NaN;
    /**
     * Column containing the allele frequency.
     */
    public String afColumn = null;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public PrsPruneOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsPruneOptions option : PrsPruneOptions.values()) {

            if (option.mandatory && !CliUtils.hasOption(aLine, option)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The summary stats file
        String filePath = CliUtils.getOptionValue(aLine, PrsPruneOptions.trainingFile);

        trainingFile = new File(filePath);

        if (!trainingFile.exists()) {

            throw new IllegalArgumentException("Training file (" + trainingFile + ") not found.");

        }

        // The ld matrix file
        filePath = CliUtils.getOptionValue(aLine, LdPruningOptions.ldMatrix);

        ldMatrixFilePath = filePath;

        // snp id column
        if (CliUtils.hasOption(aLine, PrsPruneOptions.snpId)) {

            snpIdColumn = CliUtils.getOptionValue(aLine, PrsPruneOptions.snpId);

        }

        // chr column
        if (CliUtils.hasOption(aLine, PrsPruneOptions.chrColumn)) {

            chrColumn = CliUtils.getOptionValue(aLine, PrsPruneOptions.chrColumn);

        }

        // pos column
        if (CliUtils.hasOption(aLine, PrsPruneOptions.posColumn)) {

            posColumn = CliUtils.getOptionValue(aLine, PrsPruneOptions.posColumn);

        }

        // ref column
        if (CliUtils.hasOption(aLine, PrsPruneOptions.refColumn)) {

            refColumn = CliUtils.getOptionValue(aLine, PrsPruneOptions.refColumn);

        }

        // ea column
        if (CliUtils.hasOption(aLine, PrsPruneOptions.eaColumn)) {

            eaColumn = CliUtils.getOptionValue(aLine, PrsPruneOptions.eaColumn);

        }

        // The variables
        if (CliUtils.hasOption(aLine, PrsPruneOptions.variables)) {

            String option = CliUtils.getOptionValue(aLine, PrsPruneOptions.variables);

            variables = option.split(",");
            
        }

        // The beta pattern
        if (CliUtils.hasOption(aLine, PrsPruneOptions.betaPattern)) {

            betaPattern = CliUtils.getOptionValue(aLine, PrsPruneOptions.betaPattern);

        }

        // The se pattern
        if (CliUtils.hasOption(aLine, PrsPruneOptions.sePattern)) {

            sePattern = CliUtils.getOptionValue(aLine, PrsPruneOptions.sePattern);

        }

        // The p pattern
        if (CliUtils.hasOption(aLine, PrsPruneOptions.pPattern)) {

            pPattern = CliUtils.getOptionValue(aLine, PrsPruneOptions.pPattern);

        }

        // The ld threshold
        if (CliUtils.hasOption(aLine, PrsPruneOptions.ldLocusThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsPruneOptions.ldLocusThreshold);

            try {

                ldLocusThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                throw new IllegalArgumentException("The value for LD locus threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(ldLocusThreshold) || Double.isInfinite(ldLocusThreshold)) {

                throw new IllegalArgumentException("The value for LD locus threshold (" + option + ") could not be parsed as a number.");

            }
            if (ldLocusThreshold < 0 || ldLocusThreshold > 1) {

                throw new IllegalArgumentException("The LD locus threshold (" + option + ") should be higher than 0 and lower than 1.");

            }
        }

        // The ld threshold
        if (CliUtils.hasOption(aLine, PrsPruneOptions.ldTopHitThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsPruneOptions.ldTopHitThreshold);

            try {

                ldTopHitThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                throw new IllegalArgumentException("The value for LD top hit threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(ldTopHitThreshold) || Double.isInfinite(ldTopHitThreshold)) {

                throw new IllegalArgumentException("The value for LD top hit threshold (" + option + ") could not be parsed as a number.");

            }
            if (ldTopHitThreshold < 0 || ldTopHitThreshold > 1) {

                throw new IllegalArgumentException("The LD top hit threshold (" + option + ") should be higher than 0 and lower than 1.");

            }
        }

        // The ld threshold
        if (CliUtils.hasOption(aLine, PrsPruneOptions.nSnpPerLocusThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsPruneOptions.nSnpPerLocusThreshold);

            try {

                nSnpPerLocusThreshold = Integer.parseInt(option);

            } catch (Exception e) {

                throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(nSnpPerLocusThreshold) || Double.isInfinite(nSnpPerLocusThreshold)) {

                throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") could not be parsed as a number.");

            }
            if (nSnpPerLocusThreshold < 0) {

                throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") should be higher than 0.");

            }
        }

        // The p-value threshold
        if (CliUtils.hasOption(aLine, PrsPruneOptions.pValueThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsPruneOptions.pValueThreshold);

            try {

                pValueThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("The value for p-value threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(pValueThreshold) || Double.isInfinite(pValueThreshold)) {

                throw new IllegalArgumentException("The p-value threshold (" + option + ") could not be parsed as a number.");

            }
            if (pValueThreshold <= 0 || pValueThreshold > 1) {

                throw new IllegalArgumentException("The p-value threshold (" + option + ") should be higher than 0 and lower than 1.");

            }
        }

        // The af threshold
        if (CliUtils.hasOption(aLine, PrsScoreOptions.afThreshold)) {

            String option = CliUtils.getOptionValue(aLine, PrsScoreOptions.afThreshold);

            try {

                afThreshold = Double.parseDouble(option);

            } catch (Exception e) {

                e.printStackTrace();
                throw new IllegalArgumentException("The value for allele frequency threshold (" + option + ") could not be parsed as a number.");

            }
            if (Double.isNaN(afThreshold) || Double.isInfinite(afThreshold)) {

                throw new IllegalArgumentException("The allele frequency threshold (" + option + ") could not be parsed as a number.");

            }
            if (afThreshold <= 0 || afThreshold >= 0.5) {

                throw new IllegalArgumentException("The allele frequency threshold (" + option + ") should be higher than 0 and lower than 0.5.");

            }
            
            if (!CliUtils.hasOption(aLine, PrsScoreOptions.afColumn)) {

                throw new IllegalArgumentException("No column provided for allele frequency.");
                
            }
            
            afColumn = CliUtils.getOptionValue(aLine, PrsScoreOptions.afColumn);
            
        }

        // The output file
        filePath = CliUtils.getOptionValue(aLine, PrsPruneOptions.out);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
