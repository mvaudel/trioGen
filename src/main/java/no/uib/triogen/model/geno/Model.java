package no.uib.triogen.model.geno;

import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;
import no.uib.triogen.io.IoUtils;

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
            new String[]{"β1", "β2", "β3", "β4"}
    ),
    cmf(
            "child-mother-father",
            "Regression against the number of alternative alleles of the child, mother, and father.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf_mt", "cmf_mnt", "cmf_ft", "cmf_fnt"},
            new String[]{"βc", "βm", "βf"}
    ),
    cmf_mt(
            "child-mother-father_mother-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted maternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βmt h1 + ε",
            new String[0],
            new String[]{"βc", "βm", "βf", "βmt"}
    ),
    cmf_mnt(
            "child-mother-father_mother-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and non-transmitted maternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βmnt h2 + ε",
            new String[0],
            new String[]{"βc", "βm", "βf", "βmnt"}
    ),
    cmf_ft(
            "child-mother-father_father-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and transmitted paternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βft h1 + ε",
            new String[0],
            new String[]{"βc", "βm", "βf", "βft"}
    ),
    cmf_fnt(
            "child-mother-father_father-transmitted",
            "Regression against the number of alternative alleles of the child, mother, and father, and non-transmitted pmaternal allele.",
            "y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βfnt h2 + ε",
            new String[0],
            new String[]{"βc", "βm", "βf", "βfnt"}
    ),
    cm(
            "child-mother",
            "Regression against the number of alternative alleles of the child and mother.",
            "y = βm (h1 + h2) + βc (h2 + h3) + ε",
            new String[]{"h", "cmf"},
            new String[]{"βc", "βm"}
    ),
    cf(
            "child-father",
            "Regression against the number of alternative alleles of the child and father.",
            "y = βc (h2 + h3) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf"},
            new String[]{"βc", "βf"}
    ),
    mf(
            "mother-father",
            "Regression against the number of alternative alleles of the mother and father.",
            "y = βm (h1 + h2) + βf (h3 + h4) + ε",
            new String[]{"h", "cmf"},
            new String[]{"βm", "βf"}
    ),
    c(
            "child",
            "Regression against the number of alternative alleles of the child.",
            "y = βc (h2 + h3) + ε",
            new String[]{"h", "cmf", "cm", "cf"},
            new String[]{"βc"}
    ),
    m(
            "mother",
            "Regression against the number of alternative alleles of the mother.",
            "y = βm (h2 + h3) + ε",
            new String[]{"h", "cmf", "cm", "mf"},
            new String[]{"βm"}
    ),
    f(
            "father",
            "Regression against the number of alternative alleles of the father.",
            "y = βf (h2 + h3) + ε",
            new String[]{"h", "cmf", "cf", "mf"},
            new String[]{"βf"}
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
     * @param x the x matrix
     * @param model the model
     * @param iterationI the index where to fill the matrix
     * @param h the array of hs
     * @param nAltChild the number of alternative alleles of the child
     * @param nAltMother the number of alternative alleles of the mother
     * @param nAltFather the number of alternative alleles of the father
     */
    public static void fillX(
            double[][] x,
            Model model,
            int iterationI,
            int[] h,
            int nAltChild,
            int nAltMother,
            int nAltFather
    ) {

        switch (model) {

            case h:

                x[iterationI][0] = h[0];
                x[iterationI][1] = h[1];
                x[iterationI][2] = h[2];
                x[iterationI][3] = h[3];
                return;

            case cmf:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                x[iterationI][2] = nAltFather;
                return;

            case cmf_mt:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                x[iterationI][2] = nAltFather;
                x[iterationI][3] = h[0];
                return;

            case cmf_mnt:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                x[iterationI][2] = nAltFather;
                x[iterationI][3] = h[1];
                return;

            case cmf_ft:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                x[iterationI][2] = nAltFather;
                x[iterationI][3] = h[2];
                return;

            case cmf_fnt:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                x[iterationI][2] = nAltFather;
                x[iterationI][3] = h[3];
                return;

            case cm:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltMother;
                return;

            case cf:
                x[iterationI][0] = nAltChild;
                x[iterationI][1] = nAltFather;
                return;

            case mf:
                x[iterationI][0] = nAltMother;
                x[iterationI][1] = nAltFather;
                return;

            case c:
                x[iterationI][0] = nAltChild;
                return;

            case m:
                x[iterationI][0] = nAltMother;
                return;

            case f:
                x[iterationI][0] = nAltFather;
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
            case cmf_mnt:
            case cmf_ft:
            case cmf_fnt:
                return childNotSingular && motherNotSingular && fatherNotSingular && hNotSingluar;

            case cm:
                return childNotSingular && motherNotSingular;

            case cf:
                return childNotSingular && fatherNotSingular;

            case mf:
                return motherNotSingular && fatherNotSingular;

            case c:
                return childNotSingular;

            case m:
                return motherNotSingular;

            case f:
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
     * Appends the header for the results corresponding to this model to the given string builder.
     * 
     * @param stringBuilder the string builder where to append the header
     */
    public void getHeader(
            StringBuilder stringBuilder
    ) {

        Arrays.stream(includedParentModels)
                .map(
                        parentModel -> String.join("_", name(), parentModel, "p")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.separator)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join("_", name(), betaName)
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.separator)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join("_", name(), betaName, "se")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.separator)
                                .append(label)
                );

        Arrays.stream(betaNames)
                .map(
                        betaName -> String.join("_", name(), betaName, "p")
                )
                .forEach(
                        label -> stringBuilder
                                .append(IoUtils.separator)
                                .append(label)
                );

    }
}
