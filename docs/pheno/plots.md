# TrioGen test phenotypes

Genotyped samples only, ADHD cases and ethnic outliers removed (N = 27451)


Phenotypes version V10_1.0.0-190506, standardization using [GAMLSS](https://www.gamlss.com/).


| Name | variable | Formula | Distribution | Normalization | n |
| --------- | ------- | ------------ | ------------- | - |
| Standardized Mother BMI | z_mother_bmi | `mother_bmi ~ fp(mother_age)` | `LOGNO` | `centiles.pred` Z-scores | 2192 |
| Standardized Mother Delta BMI | z_mother_delta_bmi | `mother_delta_bmi ~ cs(z_mother_bmi)` | `NO` | `centiles.pred` Z-scores | 1376 |
| Standardized Father BMI | z_father_bmi | `father_bmi ~ fp(father_age)` | `LOGNO` | `centiles.pred` Z-scores | 24545 |
| Standardized Placenta Weight | z_placenta_weight | `placenta_weight ~ fp(pregnancy_duration)` per child sex | `BCT` | `centiles.pred` Z-scores | 26700 |
| Standardized Umbilical Cord Length | z_umbilical_chord_length | `umbilical_chord_length ~ fp(pregnancy_duration)` per child sex | `BCT` | `centiles.pred` Z-scores | 26269 |
| BMI at Birth | z_bmi0 | `z_bmi0 ~ fp(pregnancy_duration)` per child sex | `LOGNO` | `centiles.pred` Z-scores | 27270 |
| BMI at 6 w | z_bmi1 | `z_bmi1 ~ fp(pregnancy_duration)` per child sex | `LOGNO` | `centiles.pred` Z-scores | 17417 |
| BMI at 3 m | z_bmi2 | `z_bmi2 ~ fp(pregnancy_duration)` per child sex | `LOGNO` | `centiles.pred` Z-scores | 23119 |
| BMI at 6 m | z_bmi3 | `z_bmi3 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 23347 |
| BMI at 8 m | z_bmi4 | `z_bmi4 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 20816 |
| BMI at 1 y | z_bmi5 | `z_bmi5 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 20773 |
| BMI at 1.5 y | z_bmi6 | `z_bmi6 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 20486 |
| BMI at 2 y | z_bmi7 | `z_bmi7 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 16106 |
| BMI at 3 y | z_bmi8 | `z_bmi8 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 16401 |
| BMI at 5 y | z_bmi9 | `z_bmi9 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 13601 |
| BMI at 7 y | z_bmi10 | `z_bmi10 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 13375 |
| BMI at 8 y | z_bmi11 | `z_bmi11 ~ pregnancy_duration` per child sex | `LOGNO` | `centiles.pred` Z-scores | 10590 |
### Mother BMI vs. Mother Age

![](docs/pheno/mother_age-mother_bmi.png)


### Standardized Mother BMI vs. Mother Age

![](docs/pheno/mother_age-z_mother_bmi.png)


### Standardized Mother BMI vs. Mother BMI

![](docs/pheno/mother_bmi-z_mother_bmi.png)


### Mother delta BMI vs. Standardized Mother BMI

![](docs/pheno/z_mother_bmi-mother_delta_bmi.png)


### Standardized Mother delta BMI vs. Standardized Mother BMI

![](docs/pheno/z_mother_bmi-z_mother_delta_bmi.png)


### Standardized Mother delta BMI vs. Mother delta BMI

![](docs/pheno/mother_delta_bmi-z_mother_delta_bmi.png)


### Father BMI vs. Father Age

![](docs/pheno/father_age-father_bmi.png)


### Standardized Father BMI vs. Father Age

![](docs/pheno/father_age-z_father_bmi.png)


### Standardized Father BMI vs. Father BMI

![](docs/pheno/father_bmi-z_father_bmi.png)


### Umbilical Cord Length vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-umbilical_chord_length.png)


### Placenta Weight vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-placenta_weight.png)


### Placenta Weight vs. Umbilical Cord Length

![](docs/pheno/umbilical_chord_length-placenta_weight.png)


### Standardized Umbilical Cord Length vs. Umbilical Cord Length

![](docs/pheno/umbilical_chord_length-z_umbilical_chord_length.png)


### Standardized Umbilical Cord Length vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_umbilical_chord_length.png)


### Standardized Placenta Weight vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_placenta_weight.png)


### Standardized Placenta Weight vs. Standardized Umbilical Cord Length

![](docs/pheno/z_umbilical_chord_length-z_placenta_weight.png)


### BMI at Birth vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi0.png)


### Standardized BMI at Birth vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi0.png)


### Standardized BMI at Birth vs. BMI at Birth

![](docs/pheno/bmi0-z_bmi0.png)


### BMI at 6 w vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi1.png)


### Standardized BMI at 6 w vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi1.png)


### Standardized BMI at 6 w vs. BMI at 6 w

![](docs/pheno/bmi1-z_bmi1.png)


### BMI at 3 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi2.png)


### Standardized BMI at 3 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi2.png)


### Standardized BMI at 3 m vs. BMI at 3 m

![](docs/pheno/bmi2-z_bmi2.png)


### BMI at 6 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi3.png)


### Standardized BMI at 6 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi3.png)


### Standardized BMI at 6 m vs. BMI at 6 m

![](docs/pheno/bmi3-z_bmi3.png)


### BMI at 8 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi4.png)


### Standardized BMI at 8 m vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi4.png)


### Standardized BMI at 8 m vs. BMI at 8 m

![](docs/pheno/bmi4-z_bmi4.png)


### BMI at 1 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi5.png)


### Standardized BMI at 1 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi5.png)


### Standardized BMI at 1 y vs. BMI at 1 y

![](docs/pheno/bmi5-z_bmi5.png)


### BMI at 1.5 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi6.png)


### Standardized BMI at 1.5 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi6.png)


### Standardized BMI at 1.5 y vs. BMI at 1.5 y

![](docs/pheno/bmi6-z_bmi6.png)


### BMI at 2 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi7.png)


### Standardized BMI at 2 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi7.png)


### Standardized BMI at 2 y vs. BMI at 2 y

![](docs/pheno/bmi7-z_bmi7.png)


### BMI at 3 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi8.png)


### Standardized BMI at 3 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi8.png)


### Standardized BMI at 3 y vs. BMI at 3 y

![](docs/pheno/bmi8-z_bmi8.png)


### BMI at 5 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi9.png)


### Standardized BMI at 5 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi9.png)


### Standardized BMI at 5 y vs. BMI at 5 y

![](docs/pheno/bmi9-z_bmi9.png)


### BMI at 7 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi10.png)


### Standardized BMI at 7 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi10.png)


### Standardized BMI at 7 y vs. BMI at 7 y

![](docs/pheno/bmi10-z_bmi10.png)


### BMI at 8 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-bmi11.png)


### Standardized BMI at 8 y vs. Pregnancy Duration

![](docs/pheno/pregnancy_duration-z_bmi11.png)


### Standardized BMI at 8 y vs. BMI at 8 y

![](docs/pheno/bmi11-z_bmi11.png)


