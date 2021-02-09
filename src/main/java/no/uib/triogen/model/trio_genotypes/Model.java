package no.uib.triogen.model.trio_genotypes;

import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.genotypes.bgen.reader.BgenVariantData;

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
     * Fills the given x matrix based on the model.
     *
     * @param x The x matrix.
     * @param model The model.
     * @param index The index where to fill the matrix.
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     * @param bgenVariantData The genotypes provider to use for this variant.
     * @param testedAlleleIndex The index of the tested allele.
     */
    public static void fillX(
            double[][] x,
            Model model,
            int index,
            String childId,
            String motherId,
            String fatherId,
            BgenVariantData bgenVariantData,
            int testedAlleleIndex
    ) {

        double[] haplotypes = bgenVariantData.getHaplotypes(
                childId,
                motherId,
                fatherId,
                index
        );
        double child = haplotypes[1] + haplotypes[2];
        double mother = haplotypes[0] + haplotypes[1];
        double father = haplotypes[2] + haplotypes[3];

        switch (model) {

            case h:

                x[index][0] = haplotypes[0];
                x[index][1] = haplotypes[1];
                x[index][2] = haplotypes[2];
                x[index][3] = haplotypes[3];
                return;

            case cmf:
                x[index][0] = child;
                x[index][1] = mother;
                x[index][2] = father;
                return;

            case cmf_mt:
                x[index][0] = child;
                x[index][1] = mother;
                x[index][2] = father;
                x[index][3] = haplotypes[1];
                return;

            case cmf_ft:
                x[index][0] = child;
                x[index][1] = mother;
                x[index][2] = father;
                x[index][3] = haplotypes[2];
                return;

            case cm:
                x[index][0] = child;
                x[index][1] = mother;
                return;

            case cm_mt:
                x[index][0] = child;
                x[index][1] = mother;
                x[index][2] = haplotypes[1];
                return;

            case cm_ft:
                x[index][0] = child;
                x[index][1] = mother;
                x[index][2] = haplotypes[2];
                return;

            case cf:
                x[index][0] = child;
                x[index][1] = father;
                return;

            case cf_mt:
                x[index][0] = child;
                x[index][1] = father;
                x[index][2] = haplotypes[1];
                return;

            case cf_ft:
                x[index][0] = child;
                x[index][1] = father;
                x[index][2] = haplotypes[2];
                return;

            case mf:
                x[index][0] = mother;
                x[index][1] = father;
                return;

            case c:
                x[index][0] = child;
                return;

            case c_mt:
                x[index][0] = child;
                x[index][1] = haplotypes[1];
                return;

            case c_ft:
                x[index][0] = child;
                x[index][1] = haplotypes[2];
                return;

            case m:
                x[index][0] = mother;
                return;

            case m_mt:
                x[index][0] = mother;
                x[index][1] = haplotypes[1];
                return;

            case f:
                x[index][0] = father;
                return;

            case f_ft:
                x[index][0] = father;
                x[index][1] = haplotypes[2];
                return;

            default:

                throw new UnsupportedOperationException("X matrix not implemented for model " + model + ".");

        }
    }

    /**
     * Returns a boolean indicating whether the matrix is likely not singular.
     *
     * @param model The model.
     * @param hMntMin The minimal number of maternal non-transmitted alleles.
     * @param hMntMax The maximal number of maternal non-transmitted alleles.
     * @param hMtMin The minimal number of maternal transmitted alleles.
     * @param hMtMax The maximal number of maternal transmitted alleles.
     * @param hFtMin The minimal number of paternal transmitted alleles.
     * @param hFtMax The maximal number of paternal transmitted alleles.
     * @param hFntMin The minimal number of paternal non-transmitted alleles.
     * @param hFntMax The maximal number of paternal non-transmitted alleles.
     *
     * @return a boolean indicating whether the matrix is likely not singular
     */
    public static boolean likelyNotSingular(
            Model model,
            double hMntMin,
            double hMntMax,
            double hMtMin,
            double hMtMax,
            double hFtMin,
            double hFtMax,
            double hFntMin,
            double hFntMax
    ) {

        switch (model) {

            case h:

                return hMntMax - hMntMin > 0.5 && hMtMax - hMtMin > 0.5 && hFtMax - hFtMin > 0.5 && hFntMax - hFntMin > 0.5;

            case cmf:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5;

            case cmf_mt:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5 && hMtMax - hMtMin > 0.5;

            case cmf_ft:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5 && hFtMax - hFtMin > 0.5;

            case cm:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5;

            case cm_mt:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && hMtMax - hMtMin > 0.5;

            case cm_ft:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && hFtMax - hFtMin > 0.5;

            case cf:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5;

            case cf_ft:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5 && hFtMax - hFtMin > 0.5;

            case cf_mt:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5 && hMtMax - hMtMin > 0.5;

            case mf:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5;

            case c:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5;

            case c_mt:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && hMtMax - hMtMin > 0.5;

            case c_ft:
                return Math.max(hMtMax, hFtMax) - Math.min(hMtMin, hFtMin) > 0.5 && hFtMax - hFtMin > 0.5;

            case m:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5;

            case m_mt:
                return Math.max(hMntMax, hMtMax) - Math.min(hMntMin, hMtMin) > 0.5 && hMtMax - hMtMin > 0.5;

            case f:
                return Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5;

            case f_ft:
                return Math.max(hFntMax, hFtMax) - Math.min(hFntMin, hFtMin) > 0.5 && hFtMax - hFtMin > 0.5;

            default:

                throw new UnsupportedOperationException("Singularity not implemented for model " + model + ".");

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
