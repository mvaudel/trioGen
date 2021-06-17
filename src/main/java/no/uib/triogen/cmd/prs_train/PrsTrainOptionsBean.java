package no.uib.triogen.cmd.prs_train;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;
import no.uib.triogen.cmd.ld_pruning.LdPruningOptions;
import no.uib.triogen.model.phenotypes.PhenotypesHandler;
import no.uib.triogen.model.trio_genotypes.Model;
import no.uib.triogen.processing.prs.PrsTrainer;
import org.apache.commons.cli.CommandLine;

/**
 * Parses and stores the command line options.
 *
 * @author Marc Vaudel
 */
public class PrsTrainOptionsBean {

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
     * Name of the effect allele column in the training file.
     */
    public Model model = Model.cmf;
    /**
     * Names of the variables.
     */
    public String[] variables = new String[]{"c", "m", "f"};
    /**
     * Pattern for the effect size column.
     */
    public String betaPattern = Model.cmf.name() + ".B" + PrsTrainer.VARIABLE_WILDCARD;
    /**
     * Pattern for the standard error column.
     */
    public String sePattern = Model.cmf.name() + ".B.se" + PrsTrainer.VARIABLE_WILDCARD;
    /**
     * Pattern for the p-value column.
     */
    public String pPattern = Model.cmf.name() + ".B.p" + PrsTrainer.VARIABLE_WILDCARD;
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
    public double ldTopHitThreshold = 0.9;
    /**
     * The minimal number of variants required for a locus.
     */
    public int nSnpPerLocusThreshold = 5;

    /**
     * Constructor. Parses the command line options and conducts minimal sanity
     * check.
     *
     * @param aLine a command line
     */
    public PrsTrainOptionsBean(
            CommandLine aLine
    ) {

        // Check that mandatory options are provided
        for (PrsTrainOptions option : PrsTrainOptions.values()) {

            if (option.mandatory && !aLine.hasOption(option.opt)) {

                throw new IllegalArgumentException("No value found for mandatory option " + option.opt + " (" + option.longOpt + ")");

            }
        }

        // The summary stats file
        String filePath = aLine.getOptionValue(PrsTrainOptions.trainingFile.opt);

        trainingFile = new File(filePath);

        if (!trainingFile.exists()) {

            throw new IllegalArgumentException("Training file (" + trainingFile.getParentFile() + ") not found.");

        }

        // The ld matrix file
        filePath = aLine.getOptionValue(LdPruningOptions.ldMatrix.opt);

        ldMatrixFilePath = filePath;

        // snp id column
        if (aLine.hasOption(PrsTrainOptions.snpId.opt)) {

            snpIdColumn = aLine.getOptionValue(PrsTrainOptions.snpId.opt);

        }

        // chr column
        if (aLine.hasOption(PrsTrainOptions.chrColumn.opt)) {

            chrColumn = aLine.getOptionValue(PrsTrainOptions.chrColumn.opt);

        }

        // pos column
        if (aLine.hasOption(PrsTrainOptions.posColumn.opt)) {

            posColumn = aLine.getOptionValue(PrsTrainOptions.posColumn.opt);

        }

        // ref column
        if (aLine.hasOption(PrsTrainOptions.refColumn.opt)) {

            refColumn = aLine.getOptionValue(PrsTrainOptions.refColumn.opt);

        }

        // ea column
        if (aLine.hasOption(PrsTrainOptions.eaColumn.opt)) {

            eaColumn = aLine.getOptionValue(PrsTrainOptions.eaColumn.opt);

        }

        // The model
        if (aLine.hasOption(PrsTrainOptions.model.opt)) {

            String option = aLine.getOptionValue(PrsTrainOptions.model.opt);

            model = Model.valueOf(option);

        }

        // The variables
        if (aLine.hasOption(PrsTrainOptions.variables.opt)) {

            String option = aLine.getOptionValue(PrsTrainOptions.variables.opt);

            variables = option.split(",");

            if (variables.length != model.betaNames.length) {

                String modelVariables = Arrays.stream(variables).collect(Collectors.joining(","));

                throw new IllegalArgumentException("Found " + variables.length + " variables (" + option + ") where " + model.betaNames.length + " (" + modelVariables + ") expected.");

            }
        }

        // The beta pattern
        if (aLine.hasOption(PrsTrainOptions.betaPattern.opt)) {

            betaPattern = aLine.getOptionValue(PrsTrainOptions.betaPattern.opt);

        }

        // The se pattern
        if (aLine.hasOption(PrsTrainOptions.sePattern.opt)) {

            sePattern = aLine.getOptionValue(PrsTrainOptions.sePattern.opt);

        }

        // The p pattern
        if (aLine.hasOption(PrsTrainOptions.pPattern.opt)) {

            pPattern = aLine.getOptionValue(PrsTrainOptions.pPattern.opt);

        }

        // The ld threshold
        if (aLine.hasOption(PrsTrainOptions.ldLocusThreshold.opt)) {

            String option = aLine.getOptionValue(PrsTrainOptions.ldLocusThreshold.opt);

            try {

                ldLocusThreshold = Double.parseDouble(option);

                if (Double.isNaN(ldLocusThreshold) || Double.isFinite(ldLocusThreshold)) {

                    throw new IllegalArgumentException("The value for LD locus threshold (" + option + ") could not be parsed as a number.");

                }
                if (ldLocusThreshold < 0 || ldLocusThreshold > 1) {

                    throw new IllegalArgumentException("The LD locus threshold (" + option + ") should be higher than 0 or lower than 1.");

                }
            } catch (Exception e) {

                throw new IllegalArgumentException("The value for LD locus threshold (" + option + ") could not be parsed as a number.");

            }
        }

        // The ld threshold
        if (aLine.hasOption(PrsTrainOptions.ldTopHitThreshold.opt)) {

            String option = aLine.getOptionValue(PrsTrainOptions.ldTopHitThreshold.opt);

            try {

                ldTopHitThreshold = Double.parseDouble(option);

                if (Double.isNaN(ldTopHitThreshold) || Double.isFinite(ldTopHitThreshold)) {

                    throw new IllegalArgumentException("The value for LD top hit threshold (" + option + ") could not be parsed as a number.");

                }
                if (ldTopHitThreshold < 0 || ldTopHitThreshold > 1) {

                    throw new IllegalArgumentException("The LD top hit threshold (" + option + ") should be higher than 0 or lower than 1.");

                }
            } catch (Exception e) {

                throw new IllegalArgumentException("The value for LD top hit threshold (" + option + ") could not be parsed as a number.");

            }
        }

        // The ld threshold
        if (aLine.hasOption(PrsTrainOptions.nSnpPerLocusThreshold.opt)) {

            String option = aLine.getOptionValue(PrsTrainOptions.nSnpPerLocusThreshold.opt);

            try {

                nSnpPerLocusThreshold = Integer.parseInt(option);

                if (Double.isNaN(nSnpPerLocusThreshold) || Double.isFinite(nSnpPerLocusThreshold)) {

                    throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") could not be parsed as a number.");

                }
                if (nSnpPerLocusThreshold < 0) {

                    throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") should be higher than 0.");

                }
            } catch (Exception e) {

                throw new IllegalArgumentException("The minimal number of SNP in locus (" + option + ") could not be parsed as a number.");

            }
        }

        // The output file
        filePath = aLine.getOptionValue(PrsTrainOptions.out.opt);

        destinationFile = new File(filePath);

        if (!destinationFile.getParentFile().exists()) {

            throw new IllegalArgumentException("Output folder (" + destinationFile.getParentFile() + ") not found.");

        }

    }
}
