# MOBA test

TrioGen v.0.3.0-beta on 27,451 full trios, ADHD cases, related, and ethnic outliers excluded. 10 PCs and genotyping batch as covariates. Reference panel AF > 5%, info score >= 0.7, singularities excluded

z_BMI “à la Chris”, see details [here](../pheno/plots.md). Note:
- not the same as the pheno tables used so far
- BMI adjusted for pregnancy duration


4 Models:

- h: `y = β1 h1 + β2 h2 + β3 h3 + β4 h4  + ε`

- cmf: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + ε`

- cmf_mt: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + βmt h1 + ε`

- cmf_ft: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + βft h1 + ε`

### Phenotypes

1. [Breastmilk Duration](breastmilk_duration.md)

1. [Formula Frequency at 6 months](formula_freq_6m.md)

1- [Pregnancy Duration](pregnancy_duration.md)

1- [BMI at birth (z_BMI0)](z_bmi0.md)

1- [BMI at 6 w (z_BMI1)](z_bmi1.md)

1- [BMI at 3 m (z_BMI2)](z_bmi2.md)

1- [BMI at 6 m (z_BMI3)](z_bmi3.md)

1- [BMI at 8 m (z_BMI4)](z_bmi4.md)

1- [BMI at 1 y (z_BMI5)](z_bmi5.md)

1- [BMI at 1.5 y (z_BMI6)](z_bmi6.md)

1- [BMI at 2 y (z_BMI7)](z_bmi7.md)

1- [BMI at 3 y (z_BMI8)](z_bmi8.md)

1- [BMI at 5 y (z_BMI9)](z_bmi9.md)

1- [BMI at 7 y (z_BMI10)](z_bmi10.md)

1- [BMI at 8 y (z_BMI11)](z_bmi11.md)

1- [BMI father (Z-score)](z_father_bmi.md)

1- [height mother (Z-score)](z_mother_height.md)

1- [Placental Weight (Z-score)](z_placenta_weight.md)

1- [Umbilical Cord Length (Z-score)](z_umbilical_chord_length.md)
