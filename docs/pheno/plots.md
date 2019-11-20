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

