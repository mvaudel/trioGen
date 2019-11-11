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
    paternal_2("h4 - h1 + h2"), 
    average_paternal("(h4 + h3 - h1 + h2)/2"), 
    child_1("h1 - h2"), 
    child_2("h3 - h4"), 
    average_child("(h1 - h2 + h3 - h4)/2");
    
    /**
     * The h nomenclature according to Chen et al. (https://doi.org/10.1101/737106).
     */
    public final String h;
    
    /**
     * Constructor.
     * 
     * @param h the h nomenclature according to Chen et al. (https://doi.org/10.1101/737106)
     */
    private TransmissionModel(
            String h
    ) {
        this.h = h;
    }

}
