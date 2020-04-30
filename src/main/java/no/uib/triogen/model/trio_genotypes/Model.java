package no.uib.triogen.model.trio_genotypes;

import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;
import no.uib.triogen.io.genotypes.GenotypesProvider;

/**
 * Enum of the linear regression models implemented.
 *
 * @author Marc Vaudel
 */
public enum Model {

    h(
            "h",
            "Regression on the h as defined by Chen et al..",
            "y = β1 h1 + β2 h2 + β3 h3 + β4 h4 + ε",
            new String[0],
            new String[]{"B1", "B2", "B3", "B4"}
    ),
    cmf(
            "child-mother-father",
            "Regression against the number of alternative alleles of the child, mother, and father.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf_mt", "cmf_mnt", "cmf_ft", "cmf_fnt"},
            new String[]{"Bc", "Bm", "Bf"}
    ),
    cmf_mt(
            "child-mother-father_mother-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted maternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βmt h1 + ε",
            new String[0],
            new String[]{"Bc", "Bm", "Bf", "Bmt"}
    ),
    cmf_ft(
            "child-mother-father_father-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted paternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βft h3 + ε",
            new String[0],
            new String[]{"Bc", "Bm", "Bf", "Bft"}
    ),
    cm(
            "child-mother",
            "Regression against the number of alternative alleles of the child and mother.",
            "y = βm (h1 + h2) + βc (h2 + h3) + ε",
            new String[]{"h", "cmf", "cm_mt"},
            new String[]{"Bc", "Bm"}
    ),
    cm_mt(
            "child-mother_mother-transmitted",
            "Regression against the number of alternative alleles of the child and mother, and transmitted maternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βmt h1 + ε",
            new String[]{"h", "cmf_mt"},
            new String[]{"Bc", "Bm", "Bmt"}
    ),
    cf(
            "child-father",
            "Regression against the number of alternative alleles of the child and father.",
            "y = βc (h2 + h3) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf", "cf_ft"},
            new String[]{"Bc", "Bf"}
    ),
    cf_ft(
            "child-father_father-transmitted",
            "Regression against the number of alternative alleles of the child and father, and transmitted paternal allele.",
            "y = βc (h2 + h3) + βf (h3 + h4) + βft h3 + ε",
            new String[]{"h", "cmf_ft"},
            new String[]{"Bc", "Bf", "Bft"}
    ),
    mf(
            "mother-father",
            "Regression against the number of alternative alleles of the mother and father.",
            "y = βm (h1 + h2) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf"},
            new String[]{"Bm", "Bf"}
    ),
    c(
            "child",
            "Regression against the number of alternative alleles of the child.",
            "y = βc (h2 + h3) + ε",
            new String[]{"h", "cmf", "cm", "cf"},
            new String[]{"Bc"}
    ),
    c_mt(
            "child_mother-transmitted",
            "Regression against the number of alternative alleles of the child and transmitted maternal allele.",
            "y = βc (h2 + h3) + βmt h1 + ε",
            new String[]{"h", "cmf_mt", "cm_mt"},
            new String[]{"Bc", "Bmt"}
    ),
    c_ft(
            "child_father-transmitted",
            "Regression against the number of alternative alleles of the child and transmitted paternal allele.",
            "y = βc (h2 + h3) + βft h3 + ε",
            new String[]{"h", "cmf_ft", "cf_ft"},
            new String[]{"Bc", "Bft"}
    ),
    m(
            "mother",
            "Regression against the number of alternative alleles of the mother.",
            "y = βm (h1 + h2) + ε",
            new String[]{"h", "cmf", "cm", "mf"},
            new String[]{"Bm"}
    ),
    m_mt(
            "mother_mother-transmitted",
            "Regression against the number of alternative alleles of the mother and transmitted maternal allele.",
            "y = βm (h1 + h2) + βmt h1 + ε",
            new String[]{"h", "cmf_mt", "cm_mt"},
            new String[]{"Bm", "Bmt"}
    ),
    f(
            "father",
            "Regression against the number of alternative alleles of the father.",
            "y = βf (h3 + h4) + ε",
            new String[]{"h", "cmf", "cf", "mf"},
            new String[]{"Bf"}
    ),
    f_ft(
            "father",
            "Regression against the number of alternative alleles of the father.",
            "y = βf (h3 + h4) + βft h3 + ε",
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

        return new String[]{"h", "cmf_mt", "cmf_ft", "cmf", "cm", "cm_mt", "c", "c_mt", "m", "m_mt"};

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
     * @param genotypesProvider The genotypes provider to use for this variant.
     * @param useDosages If true, dosages will be used when possible, hard calls
     * otherwise.
     */
    public static void fillX(
            double[][] x,
            Model model,
            int index,
            String childId,
            String motherId,
            String fatherId,
            GenotypesProvider genotypesProvider,
            boolean useDosages
    ) {
        

    }

    /**
     * Fills the given x matrix based on the model using hard calls.
     *
     * @param x The x matrix.
     * @param model The model.
     * @param index The index where to fill the matrix.
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     * @param genotypesProvider The genotypes provider to use for this variant.
     */
    public static void fillXHardCalls(
            double[][] x,
            Model model,
            int index,
            String childId,
            String motherId,
            String fatherId,
            GenotypesProvider genotypesProvider
    ) {

        short[] h = genotypesProvider.getH(
                childId,
                motherId,
                fatherId
        );
        short nAltChild = (short) (h[0] + h[2]);
        short nAltMother = (short) (h[0] + h[1]);
        short nAltFather = (short) (h[2] + h[3]);

        switch (model) {

            case h:

                x[index][0] = h[0];
                x[index][1] = h[1];
                x[index][2] = h[2];
                x[index][3] = h[3];
                return;

            case cmf:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = nAltFather;
                return;

            case cmf_mt:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = nAltFather;
                x[index][3] = h[0];
                return;

            case cmf_ft:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = nAltFather;
                x[index][3] = h[2];
                return;

            case cm:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                return;

            case cm_mt:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = h[0];
                return;

            case cf:
                x[index][0] = nAltChild;
                x[index][1] = nAltFather;
                return;

            case cf_ft:
                x[index][0] = nAltChild;
                x[index][1] = nAltFather;
                x[index][2] = h[2];
                return;

            case mf:
                x[index][0] = nAltMother;
                x[index][1] = nAltFather;
                return;

            case c:
                x[index][0] = nAltChild;
                return;

            case c_mt:
                x[index][0] = nAltChild;
                x[index][1] = h[0];
                return;

            case c_ft:
                x[index][0] = nAltChild;
                x[index][1] = h[2];
                return;

            case m:
                x[index][0] = nAltMother;
                return;

            case m_mt:
                x[index][0] = nAltMother;
                x[index][1] = h[0];
                return;

            case f:
                x[index][0] = nAltFather;
                return;

            case f_ft:
                x[index][0] = nAltFather;
                x[index][1] = h[2];
                return;

            default:

                throw new UnsupportedOperationException("X matrix not implemented for model " + model + ".");

        }
    }

    /**
     * Fills the given x matrix based on the model using dosages.
     *
     * @param x The x matrix.
     * @param model The model.
     * @param index The index where to fill the matrix.
     * @param childId The id of the child.
     * @param motherId The id of the mother.
     * @param fatherId The id of the father.
     * @param genotypesProvider The genotypes provider to use for this variant.
     */
    public static void fillXDosages(
            double[][] x,
            Model model,
            int index,
            String childId,
            String motherId,
            String fatherId,
            GenotypesProvider genotypesProvider
    ) {

        short[] h = genotypesProvider.getH(
                childId,
                motherId,
                fatherId
        );
        short nAltChild = (short) (h[0] + h[2]);
        short nAltMother = (short) (h[0] + h[1]);
        short nAltFather = (short) (h[2] + h[3]);
        
        float[] dosagesChild = genotypesProvider.getDosages(childId);
        float[] dosagesMother = genotypesProvider.getDosages(motherId);
        float[] dosagesFather = genotypesProvider.getDosages(fatherId);
        
        double dosageGenotypeChild = dosagesChild[1] + 2 * dosagesChild[2];
        double dosageGenotypeMother = dosagesMother[1] + 2 * dosagesMother[2];
        double dosageGenotypeFather = dosagesFather[1] + 2 * dosagesFather[2];

        switch (model) {

            case h:

                x[index][0] = h[0];
                x[index][1] = h[1];
                x[index][2] = h[2];
                x[index][3] = h[3];
                return;

            case cmf:
                x[index][0] = dosageGenotypeChild;
                x[index][1] = dosageGenotypeMother;
                x[index][2] = dosageGenotypeFather;
                return;

            case cmf_mt:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = nAltFather;
                x[index][3] = h[0];
                return;

            case cmf_ft:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = nAltFather;
                x[index][3] = h[2];
                return;

            case cm:
                x[index][0] = dosageGenotypeChild;
                x[index][1] = dosageGenotypeMother;
                return;

            case cm_mt:
                x[index][0] = nAltChild;
                x[index][1] = nAltMother;
                x[index][2] = h[0];
                return;

            case cf:
                x[index][0] = dosageGenotypeChild;
                x[index][1] = dosageGenotypeFather;
                return;

            case cf_ft:
                x[index][0] = nAltChild;
                x[index][1] = nAltFather;
                x[index][2] = h[2];
                return;

            case mf:
                x[index][0] = dosageGenotypeMother;
                x[index][1] = dosageGenotypeFather;
                return;

            case c:
                x[index][0] = dosageGenotypeChild;
                return;

            case c_mt:
                x[index][0] = nAltChild;
                x[index][1] = h[0];
                return;

            case c_ft:
                x[index][0] = nAltChild;
                x[index][1] = h[2];
                return;

            case m:
                x[index][0] = dosageGenotypeMother;
                return;

            case m_mt:
                x[index][0] = nAltMother;
                x[index][1] = h[0];
                return;

            case f:
                x[index][0] = dosageGenotypeFather;
                return;

            case f_ft:
                x[index][0] = nAltFather;
                x[index][1] = h[2];
                return;

            default:

                throw new UnsupportedOperationException("X matrix not implemented for model " + model + ".");

        }
    }

    /**
     * Returns a boolean indicating whether the matrix is likely not singular.
     *
     * @param model the model
     * @param hNotSingluar boolean indicating whether the h matrix is likely not
     * singular
     * @param childNotSingular boolean indicating whether the child matrix is
     * likely not singular
     * @param motherNotSingular boolean indicating whether the mother matrix is
     * likely not singular
     * @param fatherNotSingular boolean indicating whether the father matrix is
     * likely not singular
     *
     * @return a boolean indicating whether the matrix is likely not singular
     */
    public static boolean likelyNotSingular(
            Model model,
            boolean hNotSingluar,
            boolean childNotSingular,
            boolean motherNotSingular,
            boolean fatherNotSingular
    ) {

        switch (model) {

            case h:

                return hNotSingluar;

            case cmf:
                return childNotSingular && motherNotSingular && fatherNotSingular;

            case cmf_mt:
            case cmf_ft:
                return childNotSingular && motherNotSingular && fatherNotSingular && hNotSingluar;

            case cm:
            case cm_mt:
                return childNotSingular && motherNotSingular;

            case cf:
            case cf_ft:
                return childNotSingular && fatherNotSingular;

            case mf:
                return motherNotSingular && fatherNotSingular;

            case c:
            case c_mt:
            case c_ft:
                return childNotSingular;

            case m:
            case m_mt:
                return motherNotSingular;

            case f:
            case f_ft:
                return fatherNotSingular;

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
