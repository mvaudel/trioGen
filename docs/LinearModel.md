# LinearModel

This command line runs simple linear regression between phenotypes and levels of transmitted alternative alleles.

```
Y ~ h
```

For children phenotypes `Y`, TrioGen runs a linear regression with `hi`, with _i_ in {1, 2, 3, 4}, as defined by [Chen _et al._](https://doi.org/10.1101/737106).

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
-gf/--genoFormat   The genotypes file format. 0: VCF, 1: Sanger VCF.
-p/--phenoFile     The phenotypes file.
-pn/--phenoName    List of the names of the phenotypes in the phenotype file (Example: pheno1,pheno2).
-f/--fam           The trio identifiers file.
-o/--out           The file where to write the results.
```


### Additional Parameters

```
-nv/--nVariants    The number of variants to process in parallel. Default is 8.
-z/--timeOut       The number of days before timeout, default is 365..
-t/--test          If present, runs only othe first 1000 variants.
```

### Command line example

The example below runs test files provided with the tool. Please note that the command needs to be run from the folder of the repository. Note that you need to replace `Z.Y.Z` by the version number.

```
java -Xmx16G -cp bin/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.association.LinearModel -g src/main/resources/transmission/test_transmission.vcf -gf 1 -f src/main/resources/transmission/test_trio -p src/main/resources/transmission/phenos_linear_model.txt -pn pheno1,pheno2,pheno3,pheno4 -o src/test/resources/transmission/result.gz
```


### Phenotypes

Phenotypes must be provided as a tab-separated text file with phenotypes in column and samples in line, with one column named `child_id` containing children identifiers, one line per child. Phenotypes should be numeric, missing values set to `NA` or `NaN`. Please note that no normalization or standardization is conducted prior to running the linear regression.
An example of phenotype file can be found in `src/main/resources/transmission/phenos_linear_model.txt`.

The names of the columns to use for the regressions must be provided as comma-separated list, example: `pheno1,pheno2,pheno3,pheno4`. Spaces and quotes are not supported, please refrain from using spaces in column names. When multiple phenotypes are provided the regressions are conducted in parallel.


### Output

The output file contains the results of the linear regression, one line per regression, _i.e._ one per _{phenotype, variantID, h}_. The regression is conducted using the [Commons Math library](http://commons.apache.org/proper/commons-math/) and the documentation borrows information from the library documentation. For version details, please check the [pom file](https://github.com/mvaudel/trioGen/blob/master/pom.xml). 

| Column | Description |
| ------ | ----------- |
| phenotype | The name of the phenotype |
| variantID | The id of the variant |
| h | The h used for the regression |
| beta | The slope estimate |
| betaSE | Rhe [standard error of the slope estimate[(http://www.xycoon.com/standerrorb(1).htm) _s(b1)_. |
| p | The significance. |
| nH | The frequency of the h for samples with phenotype available. _E.g._ 0:23,1:2 is 23 h=0 and 2 h=1. |
| n | the number of observations that have been used for the regression. |

As detained in the [Commons Math library]](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()):
> The significance level of the slope (equiv) correlation. Specifically, the returned value is the smallest alpha such that the slope confidence interval with significance level equal to alpha does not include 0. On regression output, this is often denoted Prob(|t| > 0).
However, note that to avoid rounding of the very low p-values, the [getSignificance](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/SimpleRegression.html#getSignificance()) function is not used and instead the [regularized beta](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/special/Beta.html#regularizedBeta(double,%20double,%20double,%20double)) function is used directly.



