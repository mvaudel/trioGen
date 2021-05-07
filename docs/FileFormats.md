## Files and Formats

> **Warning: only minimal sanity check is conducted.**

- Phased genotypes in the bgen 1.3 format, one per chromosome. This can be obtained using the [_VcfToBgen_](cli/VcfToBgen.md) command line.
- Phenotypes in a tab-separated text file.
- 

### Paths

It is recommended to use absolute paths. The output folder must exist.


### Input and Output

Input files can be gzipped or not. Output files are gzipped. 

> **Warning: there is not guarantee that the order of the lines or columns in the output is the same as in the input.**

### Genotyping files

Currently, only monoallelic variants from diploid chromosomes are supported. Each variant must have a unique identifier. Genotypes of the children must be phased, with the `A|B` notation representing the allele inherited from the mother to the left (`A`) and from the father to the right (`B`). Genotyping files must be in a valid [_vcf_ format](https://www.htslib.org/doc/vcf.html). All samples should be in the same vcf files. One vcf file is processed at a time, if genotypes are spread among multiple vcf files, _e.g._ one per chromosome, use multiple command lines.

Two parsers are available and can be selected when running the command lines (see details on each command line).

1. VCF: This option uses the generic [htslib](https://www.htslib.org/) parser and can parse [all vcf files supported by htslib](https://github.com/samtools/hts-specs). The vcf file must be [compressed and indexed](https://genometoolbox.blogspot.com/2014/09/how-to-index-vcf-file.html).

2. Sanger VCF: This option uses a custom parser tested on VCF files produced by the Sanger Imputatoin service only, allowing a faster processing of the file but is not compatible with all vcf files. This parser is provided with no guarantee, and we recommend using the _Generic VCF_ parser for testing and if this parser does not work on your files. The vcf file must containt the 8 mandatory columns at the beginning. The hard calls of the genotypes of the samples must be indicated in the beginning of the field by a `0` or `1` indicating _ref_ and _alt_ alleles, respectively, and separated by a pipe (`|`) or a slash (`/`). Example: `0|1:bla:bla:bla`. Columns shoud be tab-separated. The vcf file can be gzipped or not.


### Trio file

Trios must be provided as text file. The file can be gzipped or not. The file must contain three columns with headers containing the identifiers of the children, father, and mother, in that order. Columns must be space-separated.


### Phenotypes

Phenotypes must be provided as a tab-separated text file with phenotypes in column and samples in line, with one column named `child_SentrixID` containing children identifiers, one line per child. Phenotypes should be numeric, missing values set to `NA` or `NaN`. Please note that no normalization or standardization is conducted prior to running the linear regression. The file can be gzipped.
An example of phenotype file can be found in `src/main/resources/transmission/phenos_linear_model.txt`.

The names of the columns to use for the regressions must be provided as comma-separated list, example: `pheno1,pheno2,pheno3,pheno4`. Spaces and quotes are not supported, please refrain from using spaces in column names. When multiple phenotypes are provided the regressions are conducted in parallel.


### Covariates

It is possible to include covariates in the phenotypes file. Covariates should be numeric and not contain missing or infinite values. For factors, we recommend [one-hot endcoding](https://en.wikipedia.org/wiki/One-hot). These covariates are typically PCs or batches, we recommend to account for phenotypic covariates like age when normalizing the phenotypes, see [gamlss](https://www.gamlss.com/) for examples.

If covariates are provided, the [OLS implementation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/stat/regression/OLSMultipleLinearRegression.html) of the [Commons Math library](http://commons.apache.org/proper/commons-math/) is used to run a linear regression between the phenotypes and the matrix of covariates.

```
y = ßi xi + e                                                                             (covariates)
```

where y is the phenotypes, xi the ith covariates, and e the residuals. The residuals are then used as new phenotypes.

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

