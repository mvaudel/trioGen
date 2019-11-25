# LinearModel

This command line runs linear regression between phenotypes and levels of transmitted alternative alleles. As detailed below, different models are available. They are based on hard-called genotypes, and need children genotypes to be phased. For each model, the estimated regression coefficients are returned, and models are compared to infer the relevance of including mother, father, and child in the regression.

### Regression against the number of transmitted alternative alleles

In the `h` model, the regression is conducted against the number of transmitted and non-transmitted alternative alleles, h, as defined by [Chen _et al._](https://doi.org/10.1101/737106).

```
y = β1 h1 + β2 h2 + β3 h3 + β4 h4 + ε                                                     (h)
```

Where y represents the phenotypes, h1 the number of maternal alternative alleles transmitted to the child, h2 the number of maternal alternative alleles non-transmitted to the child, h3 the number of paternal alternative alleles transmitted to the child, and h4 the number of paternal alternative alleles non-transmitted to the child.

### Regression against the number of alternative alleles for the child, mother, and father 

In the `h` model above: (1) `β1 - β2` and `β3 - β4` represent the _child genetic effect_, `βc`, of the alternative alleles transmitted by the mother and father, respectively. Under the assumption that there is no difference in the child genetic effect between the allele transmitted by the mother or the one transmitted by the father, we have `βc = β1 - β2 = β3 - β4`; (2) `β2` and `β1 - βc` represent the _mother genetic effect_, `βm`, of the alternative alleles non-transmitted and transmitted by the mother, respectively. Under the assumption that there is no difference in the mother genetic effect between the transmitted and non-transmitted allele, we have `βm = β2 = β1 - βc`; and (3) `β4` and `β3 - βc` represent the _father genetic effect_, `βf`, of the alternative alleles non-transmitted and transmitted by the father, respectively. Under the assumption that there is no difference in the father genetic effect between the transmitted and non-transmitted allele, we have `βf = β4 = β3 - βc`.

Under these assumptions, the `h` model above can be written as child-mother-father `cmf` as below:

```
y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + ε                                        (cmf)
```

### Transmitted or non-transmitted maternal or paternal association

In the presence of an association with the transmitted or non-transmitted alleles specifically, `β2 ≠ β1 - βc` and `β4 ≠ β3 - βc`, for maternal and paternal alleles, respectively. This can be modeled by including `βmt`, `βft`, the regression coefficients corresponding to the maternal and paternal transmitted alleles, respectively, and their non-transmitted counterparts, `βmnt` and `βfnt`, respectively.

| Name | Variable | Definition |
| ---- | -------- | ---------- |
| Mother transmitted allele | `βmt` | `βmt = β1 - βc - βm` |
| Father transmitted allele | `βft` | `βft = β3 - βc - βf` |
| Mother non-transmitted allele | `βmnt` | `βmt = β2 - βm` |
| Father non-transmitted allele | `βfnt` | `βft = β4 - βf` |

From the previous models, we hence derive the following models:

| Name | Variable | Hypothesis |
| ---- | -------- | ---------- |
| child-mother-father_mother-transmitted | `cmf_mt` | `β1 - βc - βm ≠ 0` |
| child-mother-father_mother-non-transmitted | `cmf_mnt` | `β2 - βm ≠ 0` |
| child-mother-father_father-transmitted | `cmf_ft` | `β3 - βc - βf ≠ 0` |
| child-mother-father_father-non-transmitted | `cmf_fnt` | `β4 - βf ≠ 0` |

```
y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βmt h1 + ε                               (cmf_mt)
```

```
y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βmnt h2 + ε                              (cmf_mnt)
```

```
y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βft h3 + ε                               (cmf_ft)
```

```
y = βm (h1 + h2) + βc (h2 + h3) + βf (h3 + h4) + βfnt h4 + ε                              (cmf_fnt)
```


### Regression against the number of alternative alleles for the child, mother, or father 

From the cmf model, we can derive the following models:

| Name | Variable | Hypothesis |
| ---- | -------- | ---------- |
| child-mother | `cm` | `βf = 0` |
| child-father | `cf` | `βm = 0` |
| mother-father | `mf` | `βc = 0` |
| child | `c` | `βm = 0` and `βf = 0` |
| mother | `m` | `βc = 0` and `βf = 0` |
| father | `f` | `βc = 0` and `βm = 0` |

```
y = βm (h1 + h2) + βc (h2 + h3) + ε                                                       (cm)
```

```
y = βc (h2 + h3) + βf (h3 + h4) + ε                                                       (cf)
```

```
y = βm (h1 + h2) + βf (h3 + h4) + ε                                                       (mf)
```

```
y = βc (h2 + h3) + ε                                                                      (c)
```

```
y = βm (h1 + h2) + ε                                                                      (m)
```

```
y = βf (h3 + h4) + ε                                                                      (f)
```


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.association.LinearModel [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


### Standard parameters

```
-h/--help          Display help text
-v/--version       Display version
```


### Mandatory Parameters

```
-g/--geno          The genotypes file.
-p/--phenoFile     The phenotypes file.
-pn/--phenoName    List of the names of the phenotypes in the phenotype file (Example: pheno1,pheno2).
-f/--fam           The trio identifiers file.
-o/--out           The file where to write the results.
```


### Additional Parameters

```
-cv/--covariate    List of the names of the covariates columns in the phenotypes file. Example: pc1,pc2,pc3,pc4,pc5,pc6.
-m/--model         List of the names of the models to use. Default: h,cmf.
-v/variantId       File listing the variants to include in the analysis. Example: rs123,rs456. Default: iterate though all variants.
-gf/--genoFormat   The genotypes file format. 0: VCF, 1: Sanger VCF. Default is VCF.
-nv/--nVariants    The number of variants to process in parallel. Default is 8.
-z/--timeOut       The number of days before timeout, default is 365.
-t/--test          If present, runs only othe first 1000 variants.
```

### Command line example

The example below runs simulated test files. Please note that the command needs to be run from the folder of the repository and that you need to replace `X.Y.Z` by the version number.

```
java -Xmx16G -cp bin/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.association.LinearModel -g src/main/resources/transmission/test_transmission.vcf -gf 1 -f src/main/resources/transmission/test_trio -p src/main/resources/transmission/phenos_linear_model.txt -pn pheno1,pheno2,pheno3,pheno4 -o src/test/resources/transmission/result.gz
```


### Phenotypes

Phenotypes must be provided as a tab-separated text file with phenotypes in column and samples in line, with one column named `child_SentrixID` containing children identifiers, one line per child. Phenotypes should be numeric, missing values set to `NA` or `NaN`. Please note that no normalization or standardization is conducted prior to running the linear regression. The file can be gzipped.
An example of phenotype file can be found in `src/main/resources/transmission/phenos_linear_model.txt`.

The names of the columns to use for the regressions must be provided as comma-separated list, example: `pheno1,pheno2,pheno3,pheno4`. Spaces and quotes are not supported, please refrain from using spaces in column names. When multiple phenotypes are provided the regressions are conducted in parallel.


### Output

The output file contains the results of the linear regression, one line per phenotype per variant. The regression is conducted using the [OLS implementation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/OLSMultipleLinearRegression.html) of the [Commons Math library](http://commons.apache.org/proper/commons-math/) and the documentation borrows information from the library documentation. For version details, please check the [pom file](https://github.com/mvaudel/trioGen/blob/master/pom.xml). 

- Each line starts with information on the phenotype, variant, and allele distribution among the trios included in the regression.

| Column | Description |
| ------ | ----------- |
| `phenotype` | The name of the phenotype. |
| `variantID` | The id of the variant. |
| `n` | The number of samples included in the regression. |
| `nAlt` | The distribution of alternative alleles among the trios included in the regression. |
| `nH` | The distribution of alternative allele transmission, h, as defined by [Chen _et al._](https://doi.org/10.1101/737106) among the trios included in the regression. |

- Subsequently, for each model, estimates and summary statistics are listed.

| Column Name Scheme | Example | Description |
| ------------------ | ------- | ----------- |
| `model1_model2_p` | `mf_cmf_p` | Estimate of whether model2 provides significantly better fit to the data than model1 using an F-test. |
| `model_variable_beta` | `mf_βm` | Estimate of the slope of the variable in the model. |
| `model_variable_beta_se` | `h_β2_se` | [Standard error of the slope estimate](http://www.xycoon.com/standerrorb(1).htm). |
| `model_variable_beta_p` | `cmf_βc_p` | Significance level of the slope (equiv) correlation. |

- F-test

F is estimated using the residual sum of squares of the models and their degrees of freedom (number of slopes + intercept). The cumulative distribution function (CDF) is obtained from the [F-distribution](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/distribution/FDistribution.html) implementation of the [Commons Math library](http://commons.apache.org/proper/commons-math/).

- Significance level of the slope

As detained in the [Commons Math library](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()), the significance p corresponds to the significance level of the slope (equiv) correlation. Specifically, the returned value is the smallest alpha such that the slope confidence interval with significance level equal to alpha does not include 0. On regression output, this is often denoted Prob(|t| > 0). Note that to avoid rounding of the very low p-values, the [getSignificance](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()) function is not used and instead the [regularized beta](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/special/Beta.html#regularizedBeta(double,%20double,%20double,%20double)) function is used directly.


### Target variants

If you want to process variants specifically, please provide them as text file using the `-v/variantID` command line argument. The first lines of the text file starting with '#' will be ignored. The file must contain a single-line header that is not starting with '#'. The file must be tab-separated and containt in the four first columns: (1) the id of the variant as present in the vcf file; (2) the chromosome name; (3) the position on the chromosome where to start looking for the variant; and (4) the position on the chromosome where to stop looking for the variant.

#### Example:
```
# Variant target example
id	chromosome	start	end
rs123	1	123456	123456
rs456	2	456789	456798
rs789	3	789789	7897890
```



