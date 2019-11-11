package no.uib.triogen.model.geno;

/**
 *
 * @author Marc Vaudel
 */
public enum TransmissionModel {

    maternal_transmitted("h1"),
    maternal_non_transmitted("h2"),
    paternal_transmitted("h3"),
    paternal_non_transmitted("h4"),
    average_maternal_no_father("(h1 - h3 + h2)/2"),
    average_child_no_father("(h1 - h2 + h3)/2"),
    maternal_1("h2"),
    maternal_2("h1 - h3 + h4"),
    average_maternal("(h2 + h1 - h3 + h4)/2"),
    paternal_1("h4"),
    paternal_2("h3 - h1 + h2"),
    average_paternal("(h4 + h3 - h1 + h2)/2"),
    child_1("h1 - h2"),
    child_2("h3 - h4"),
    average_child("(h1 - h2 + h3 - h4)/2");

    /**
     * The h nomenclature according to Chen et al.
     * (https://doi.org/10.1101/737106).
     */
    public final String h;

    /**
     * Constructor.
     *
     * @param h the h nomenclature according to Chen et al.
     * (https://doi.org/10.1101/737106)
     */
    private TransmissionModel(
            String h
    ) {
        this.h = h;
    }

    /**
     * Returns the value to use in the regression.
     * 
     * @param hs the estimates of the hs
     * 
     * @return the value to use in the regression
     */
    public double getRegressionValue(double[] hs) {

        switch (this) {
            case maternal_transmitted:
                return hs[0];

            case maternal_non_transmitted:
                return hs[1];

            case paternal_transmitted:
                return hs[2];

            case paternal_non_transmitted:
                return hs[3];

            case average_maternal_no_father:
                return (hs[0] - hs[2] + hs[1]) / 2;

            case average_child_no_father:
                return (hs[0] - hs[1] + hs[2]) / 2;

            case maternal_1:
                return hs[1];

            case maternal_2:
                return hs[0] - hs[2] + hs[3];

            case average_maternal:
                return (hs[1] + hs[0] - hs[2] + hs[3]) / 2;

            case paternal_1:
                return hs[3];

            case paternal_2:
                return hs[2] - hs[0] + hs[1];

            case average_paternal:
                return (hs[3] + hs[2] - hs[0] + hs[1]) / 2;

            case child_1:
                return hs[0] - hs[1];

            case child_2:
                return hs[3] - hs[4];

            case average_child:
                return (hs[0] - hs[1] + hs[3] - hs[4]) / 2;

            default:
                throw new UnsupportedOperationException("Regression value for transmission model " + this + " not implemented.");
        }
    }
}
