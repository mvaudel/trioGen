## LinearModel

This command line runs linear regression models between phenotypes and {child, mother, father} trio genotypes. 

### General considerations

Regression can be conducted using different models, as detailed [here](../Models.md).

Files need to be formatted as detailed [here](../FileFormats.md).

The tested allele is selected as detailed [here](../AlleleFrequency.md).

Mendelian errors are handled as detailed [here](../MendelianErrors.md).

Covariates are handled as detailed [here](../Covariates.md).

Allele transmission is handled as detailed [here](../Transmission.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.association.LinearModel [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-g/--geno                 The genotypes file.
-c/--chromosome           The name of the chromosome (1-22, X/23, or Y/24 - please keep the naming of sex chromosome consistent across all files)
-p/--phenoFile            The phenotypes and covariates file.
-pn/--phenoName           List of the names of the phenotypes in the phenotype file (Example: pheno1,pheno2).
-f/--fam                  The trio identifiers file.
-o/--out                  The file where to write the results.
```


#### Additional Parameters

```
-cg/--covariate_general   List of the names of the covariates to use for all phenotypes. Example: pc1,pc2,pc3,pc4,pc5,pc6.
-cs/--covariate_specific  File containing the names of the covariates to use for a specific phenotype in json format.
-m/--model                List of the names of the models to use. Default: h,cmf_mt,cmf.
-v/--variantId            File listing the variants to include in the analysis. Default: process all variants in the genotypes file.
-d/--dist                 The maximum distance in bp to consider around a variant. Default: 500000.
-af/--afThreshold         Allele frequency threshold. 0.005 excludes all alleles of variants with frequency < 0.5% or > 99.5%. Default: 0.005.
-id/--childId             The name of the column containing the child id. Default: child_SentrixID.
-nv/--nVariants           The number of variants to process in parallel. Default is 8.
-x0/--x0                  If present the association results will only be reported when multiple values of x are available for the regression.
-z/--timeOut              The number of days before timeout, default is 365.
-vl/--variantLog          If present, writes a log for every variant next to the results file.
```


### Processing

The regression is conducted using the [OLS implementation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/OLSMultipleLinearRegression.html) of the [Commons Math library](http://commons.apache.org/proper/commons-math/) and the documentation borrows information from the library documentation. For version details, please check the [pom file](https://github.com/mvaudel/trioGen/blob/master/pom.xml). 

Effect size, [standard error](http://www.xycoon.com/standerrorb(1).htm), and significance for each variable of the model. The performance of the model is then compared to a simple intercept, and to all the models used that use the same variables but with more degrees of freedom using an F-test.

- F-test

F is estimated using the residual sum of squares of the models and their degrees of freedom (number of variables + intercept). The cumulative distribution function (CDF) is obtained from the [F-distribution](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/distribution/FDistribution.html) implementation of the [Commons Math library](http://commons.apache.org/proper/commons-math/).

- Significance level of the slope

As detained in the [Commons Math library](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()), the significance p corresponds to the significance level of the slope (equiv) correlation. Specifically, the returned value is the smallest alpha such that the slope confidence interval with significance level equal to alpha does not include 0. On regression output, this is often denoted Prob(|t| > 0). Note that to avoid rounding of the very low p-values, the [getSignificance](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()) function is not used and instead the [regularized beta](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/special/Beta.html#regularizedBeta(double,%20double,%20double,%20double)) function is used directly.


### Output

The output file contains the results of the linear regression, one line per phenotype per variant. The output is gz compressed, with gz blocks per line, and a file containing the index of each line is written next to it '*.index.gz' file. The files can be read by standard libraries, however, they tend to be large. If you need to extract specific lines or columns, please consider the [_Extract_](Extract.md) command line.

- Each line starts with information on the phenotype, variant, and allele distribution among the trios included in the regression.

| Column | Description |
| ------ | ----------- |
| `phenotype` | The name of the phenotype. |
| `contig` | The contig/chromosome containing the variant. |
| `position` | The position of the variant. |
| `variantId` | The variant identifier in the bgen file. |
| `rsId` | The rsId in the bgen file. |
| `testedAllele` | The tested allele. |
| `otherAllele` | The other allele. |
| `n` | Number of trios included in the regression. |
| `nAlt` | Number of tested alleles in the children, mothers, and fathers. |
| `nH` | Number of tested alleles in the haplotypes. |
| `mendelianError` | Estimate of the share of mendelian errors. |

- Subsequently, for each model, estimates and summary statistics are listed.

| Column Name Scheme | Example | Description |
| ------------------ | ------- | ----------- |
| `model.intercept.p` | `cmf_mt.p` | Estimate of how well the model performed compared to a simple intercept (F-test p-value). |
| `model.otherModel.p` | `cmf.h.p` | Estimate of how well the model performed compared to another model with more degrees of freedom (F-test p-value). |
| `model.variable` | `cmf.Bm` | Estimate of the effect size for the given variable. |
| `model.variable.se` | `cmf.Bm.se` | Standard error of the effect size estimate for the given variable. |
| `model.variable.p` | `cmf.Bm.p` | Significance level for the given variable. |

