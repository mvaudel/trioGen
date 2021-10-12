package no.uib.triogen.model.trio_genotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.genotypes.bgen.variant_data.BgenVariantTrioData;

/**
 * Enum of the linear regression models implemented.
 *
 * @author Marc Vaudel
 */
public enum Model {

    h(
            "h",
            "Regression on the haplotyes.",
            "y = βmnt hmnt + βmt hmt + βft hft + βfnt hfnt + ε",
            new String[0],
            new String[]{"Bmnt", "Bmt", "Bft", "Bfnt"}
    ),
    cmf(
            "child-mother-father",
            "Regression against the number of alternative alleles of the child, mother, and father.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + βf (hft + hfnt) + ε",
            new String[]{"h", "cmf_mt", "cmf_mnt", "cmf_ft", "cmf_fnt"},
            new String[]{"Bc", "Bm", "Bf"}
    ),
    cmf_mt(
            "child-mother-father_mother-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted maternal allele.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + βf (hft + hfnt) + βmt hmt + ε",
            new String[0],
            new String[]{"Bc", "Bm", "Bf", "Bmt"}
    ),
    cmf_ft(
            "child-mother-father_father-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted paternal allele.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + βf (hft + hfnt) + βft hft + ε",
            new String[0],
            new String[]{"Bc", "Bm", "Bf", "Bft"}
    ),
    cm(
            "child-mother",
            "Regression against the number of alternative alleles of the child and mother.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + ε",
            new String[]{"h", "cmf", "cm_mt", "cm_ft"},
            new String[]{"Bc", "Bm"}
    ),
    cm_mt(
            "child-mother_mother-transmitted",
            "Regression against the number of alternative alleles of the child and mother, and transmitted maternal allele.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + βmt hmt + ε",
            new String[]{"h", "cmf_mt"},
            new String[]{"Bc", "Bm", "Bmt"}
    ),
    cm_ft(
            "child-mother_mother-transmitted",
            "Regression against the number of alternative alleles of the child and mother, and transmitted maternal allele.",
            "y = βm (hmnt + hmt) + βc (hmt + hft) + βft hft + ε",
            new String[]{"h", "cmf_ft"},
            new String[]{"Bc", "Bm", "Bft"}
    ),
    cf(
            "child-father",
            "Regression against the number of alternative alleles of the child and father.",
            "y = βc (hmt + hft) + βf (hft + hfnt) + ε",
            new String[]{"h", "cmf", "cf_mt", "cf_ft"},
            new String[]{"Bc", "Bf"}
    ),
    cf_mt(
            "child-father_father-transmitted",
            "Regression against the number of alternative alleles of the child and father, and transmitted paternal allele.",
            "y = βc (hmt + hft) + βf (hfnt + hft) + βmt hmt + ε",
            new String[]{"h", "cmf_mt"},
            new String[]{"Bc", "Bf", "Bmt"}
    ),
    cf_ft(
            "child-father_father-transmitted",
            "Regression against the number of alternative alleles of the child and father, and transmitted paternal allele.",
            "y = βc (hmt + hft) + βf (hft + hfnt) + βft hft + ε",
            new String[]{"h", "cmf_ft"},
            new String[]{"Bc", "Bf", "Bft"}
    ),
    mf(
            "mother-father",
            "Regression against the number of alternative alleles of the mother and father.",
            "y = βm (hmnt + hmt) + βf (hft + hfnt) + ε",
            new String[]{"h", "cmf"},
            new String[]{"Bm", "Bf"}
    ),
    c(
            "child",
            "Regression against the number of alternative alleles of the child.",
            "y = βc (hmt + hft) + ε",
            new String[]{"h", "cmf", "cm", "cf"},
            new String[]{"Bc"}
    ),
    c_mt(
            "child_mother-transmitted",
            "Regression against the number of alternative alleles of the child and transmitted maternal allele.",
            "y = βc (hmt + hft) + βmt hmt + ε",
            new String[]{"h", "cmf_mt", "cm_mt"},
            new String[]{"Bc", "Bmt"}
    ),
    c_ft(
            "child_father-transmitted",
            "Regression against the number of alternative alleles of the child and transmitted paternal allele.",
            "y = βc (hmt + hft) + βft hft + ε",
            new String[]{"h", "cmf_ft", "cf_ft"},
            new String[]{"Bc", "Bft"}
    ),
    m(
            "mother",
            "Regression against the number of alternative alleles of the mother.",
            "y = βm (hmnt + hmt) + ε",
            new String[]{"h", "cmf", "cm", "mf"},
            new String[]{"Bm"}
    ),
    m_mt(
            "mother_mother-transmitted",
            "Regression against the number of alternative alleles of the mother and transmitted maternal allele.",
            "y = βm (hmnt + hmt) + βmt hmt + ε",
            new String[]{"h", "cmf_mt", "cm_mt"},
            new String[]{"Bm", "Bmt"}
    ),
    f(
            "father",
            "Regression against the number of alternative alleles of the father.",
            "y = βf (hft + hfnt) + ε",
            new String[]{"h", "cmf", "cf", "mf"},
            new String[]{"Bf"}
    ),
    f_ft(
            "father",
            "Regression against the number of alternative alleles of the father.",
            "y = βf (hft + hfnt) + βft hft + ε",
            new String[]{"h", "cmf_ft", "cf_ft"},
            new String[]{"Bf", "Bft"}
    );

    /**
     * The full name of the model.
     */
    public final String fullName;
    /**
     * The description of the model.
     */
    public final String description;
    /**
     * The equation of the model.
     */
    public final String equation;
    /**
     * The models this one is nested in.
     */
    public final String[] parentModels;
    /**
     * The names of the betas to estimate.
     */
    public final String[] betaNames;
    /**
     * Cache for the parent models to include in a given analysis. Needs to be
     * set using setParentModels.
     */
    public String[] includedParentModels = null;

    /**
     * Constructor.
     *
     * @param index the index of the option
     * @param fullName the full name of the option
     * @param description the description of the option
     * @param parentModels the models this one is nested in
     * @param betaNames the names of the betas to estimate
     */
    private Model(
            String fullName,
            String description,
            String equation,
            String[] parentModels,
            String[] betaNames
    ) {

        this.fullName = fullName;
        this.description = description;
        this.equation = equation;
        this.parentModels = parentModels;
        this.betaNames = betaNames;

    }

    /**
     * Returns all options as a String.
     *
     * @return all options as a String
     */
    public static String getCommandLineOptions() {

        return Arrays.stream(values())
                .map(
                        model -> model.name()
                )
                .collect(
                        Collectors.joining(",")
                );

    }

    /**
     * Returns the models to run by default.
     *
     * @return the models to run by default
     */
    public static String[] getDefaultOption() {

        return new String[]{"h", "cmf_mt", "cmf_ft", "cmf"};

    }

    /**
     * Indicates whether there is enough data to compute the given model.
     *
     * @param model The model.
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     * @param bgenVariantData The genotypes provider to use for this variant.
     *
     * @return A boolean indicating whether there is enough data to compute the
     * given model.
     */
    public static boolean hasData(
            Model model,
            String childId,
            String motherId,
            String fatherId,
            BgenVariantTrioData bgenVariantData
    ) {

        switch (model) {

            case h:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId) && bgenVariantData.contains(fatherId);

            case cmf:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId) && bgenVariantData.contains(fatherId);

            case cmf_mt:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId) && bgenVariantData.contains(fatherId);

            case cmf_ft:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId) && bgenVariantData.contains(fatherId);

            case cm:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId);

            case cm_mt:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId);

            case cm_ft:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId);

            case cf:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(fatherId);

            case cf_mt:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(fatherId);

            case cf_ft:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(fatherId);

            case mf:
                return bgenVariantData.contains(motherId) && bgenVariantData.contains(fatherId);

            case c:
                return bgenVariantData.contains(childId);

            case c_mt:
                return bgenVariantData.contains(childId);

            case c_ft:
                return bgenVariantData.contains(childId);

            case m:
                return bgenVariantData.contains(motherId);

            case m_mt:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(motherId);

            case f:
                return bgenVariantData.contains(fatherId);

            case f_ft:
                return bgenVariantData.contains(childId) && bgenVariantData.contains(fatherId);

            default:

                throw new UnsupportedOperationException("X matrix not implemented for model " + model + ".");

        }
    }

    /**
     * Returns the x matrix to use in the model.
     *
     * @param model The model.
     * @param i The row index in the matrices.
     * @param j The column index in the model.
     * @param haplotypeX The matrix of haplotypes.
     * @param childX The matrix of child genotypes.
     * @param motherX The matrix of mother genotypes.
     * @param fatherX The matrix of father genotypes.
     * 
     * @return the matrix to use in the regression.
     */
    public static double getXValueAt(
            Model model,
            int i,
            int j,
            double[][] haplotypeX,
            double[][] childX,
            double[][] motherX,
            double[][] fatherX
    ) {

        switch (model) {

            case h:
                
                return haplotypeX[i][j];

            case cmf:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    case 2: return fatherX[i][0];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cmf_mt:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    case 2: return fatherX[i][0];
                    case 3: return haplotypeX[i][1];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cmf_ft:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    case 2: return fatherX[i][0];
                    case 3: return haplotypeX[i][2];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cm:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cm_mt:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    case 2: return haplotypeX[i][1];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cm_ft:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return motherX[i][0];
                    case 2: return haplotypeX[i][2];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cf:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return fatherX[i][0];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cf_mt:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return fatherX[i][0];
                    case 2: return haplotypeX[i][1];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case cf_ft:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return fatherX[i][0];
                    case 2: return haplotypeX[i][2];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                   
                }

            case mf:
                
                switch (j) {
                    
                    case 0: return motherX[i][0];
                    case 1: return fatherX[i][0];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case c:
                
                return childX[i][0];

            case c_mt:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return haplotypeX[i][1];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case c_ft:
                
                switch (j) {
                    
                    case 0: return childX[i][0];
                    case 1: return haplotypeX[i][2];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case m:
                
                return motherX[i][0];

            case m_mt:
                
                switch (j) {
                    
                    case 0: return motherX[i][0];
                    case 1: return haplotypeX[i][1];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            case f:
                
                return fatherX[i][0];

            case f_ft:
                
                switch (j) {
                    
                    case 0: return fatherX[i][0];
                    case 1: return haplotypeX[i][2];
                    default: throw new IllegalArgumentException("Only " + model.betaNames.length + " columns in model " + model + ".");
                    
                }

            default:

                throw new UnsupportedOperationException("X matrix not implemented for model " + model + ".");

        }
    }

    /**
     * Sets the parent models to consider for this model based on the given
     * included models.
     *
     * @param includedModels the models included
     */
    public void setParentModels(
            HashSet<String> includedModels
    ) {

        includedParentModels = Arrays.stream(parentModels)
                .filter(
                        parentModel -> includedModels.contains(parentModel)
                )
                .toArray(String[]::new);

    }

    /**
     * Appends the header for the results corresponding to this model to the
     * given string builder.
     *
     * @param stringBuilder the string builder where to append the header
     */
    public void getHeader(
            StringBuilder stringBuilder
    ) {

        stringBuilder
                .append(IoUtils.SEPARATOR)
                .append(
                        String.join(".", name(), "variance_explained")
                )
                .append(IoUtils.SEPARATOR)
                .append(
                        String.join(".", name(), "intercept", "p")
                );

        Arrays.stream(includedParentModels)
                .map(
                        parentModel -> String.join(".", name(), parentModel, "p")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join(".", name(), betaName)
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join(".", name(), betaName, "se")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join(".", name(), betaName, "p")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.SEPARATOR)
                                .append(label)
                );

    }
}
