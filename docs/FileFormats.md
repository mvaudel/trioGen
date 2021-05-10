## Files and Formats

> **Warning: only minimal sanity check is conducted.**

> **Warning: for multithreaded operations, please expect some stochasticity in the order of the lines in the output, and a different order between runs.**

### Paths

It is recommended to use absolute paths. The output folder must exist.


### Compression

Input files can be gzipped or not. Output files are gzipped. 


### Genotyping files

Genotypes files must be in the bgen 1.3 format, one per chromosome, with the haplotypes phased and compressed using zstd. All samples must be in the same file. This can be obtained using the [_VcfToBgen_](cli/VcfToBgen.md) command line.

Each variant must have a unique identifier, different variants can have the same rsid.

Indels are currently not supported. If this is something important for your research, please [open an issue](https://github.com/mvaudel/trioGen/issues).

Multi-allelic variants and sex chromosomes are supported. See documentation on [_AlleleFrequency_](AlleleFrequency.md) and [_Transmission_](Transmission.md) for more details.

For haplotype analyses, the genotypes of the children must be phased. The quality of the phasing can be monitored using the estimated prevalence of Mendelian errors, as detailed [here](MendelianErrors.md).


### Trio file

The identifiers of children, mothers, and fathers trios must be provided as text file. The file can be gzipped or not. The file must contain three columns with headers containing the identifiers of the children, father, and mother, in that order. Columns must be space-separated.

When two trios share related individuals, we recommend to keep only one, and hence only work with unrelated trios.


### Covariates and Phenotypes

Covariates and phenotypes must be provided as a single tab-separated text file with phenotypes in column and samples in line, with one column containing children identifiers, one line per trio. Phenotypes must be numeric, missing values set to `NA` or `NaN`. The file can be gzipped. Please do not use spaces, quotes, or special characters in column names or identifiers.


### Specific covariates

When processing multiple phenotypes in parallel, it is possible to provide phenotype-specific covariates. These must be provided in a json format. For each phenotype name, please provide an array of covariates under the label 'controlVariables'. Other fields will be ignored.

#### Example:
```
{
	"pheno1":
	{
		"controlVariables":["covariate1","covariate2"]
	},
	"pheno2":
	{
		"controlVariables":["covariate3","covariate4"]
	}
}
```


### Trio LD Matrices

Trio LD matrices '.tld' must be provided per chromosome and generated using the [LdMatrix](cli/LdMatrix.md) command line. Please note that generating these files can take several days, but needs to be done only once. In the MoBa cohort, these files range from 200 MB to 1.2 GB in size.


### Target variants

When processing specific variants, they must be provided as a tab-separated text file. The file can be gzipped or not. The first lines of the text file starting with '#' will be ignored. The file must contain a single-line header that is not starting with '#'. The file must be tab-separated and contain the five columns: (1) the id of the variant as present in the bgen file, please refrain from using rsids if not unique; (2) the chromosome name; (3) the position on the chromosome where to start looking for the variant; and (4) the position on the chromosome where to stop looking for the variant; and (5) a column for later reference. More columns can be added and will be ignored. The order of the lines has no importance.

This file can be generated from a list of identifiers using the [VariantFile](cli/VariantFile.md) command.

#### Example:
```
# Variant target example
id	chromosome	start	end	source
rs123	1	123456	123456	pmid123
rs456	2	456789	456798	BMI
rs789	3	789789	7897890	T2D
```

> **Warning: the name of the sex chromosome must be the same as in the bgen file.**


### Simple risk scores

The weights used to compute a simple risk score must be provided as a tab-separated text file. The file can be gzipped or not. The first lines of the text file starting with '#' will be ignored. The file must contain a single-line header that is not starting with '#'. The file must be tab-separated and contain the five columns: (1) the id of the variant as present in the gben file, please refrain from using rsids if not unique; (2) the chromosome name; (3) the position on the chromosome ; and (4) the effect allele as a String; and (5) the weight as a number. More columns can be added and will be ignored. The order of the lines has no importance.

#### Example:
```
# A simple score
id	chromosome	bp	ea	weight
rs123	1	123456	A	0.1
rs456	2	456789	CG	2
rs789	3	789789	T	0.01
```

> **Warning: the name of the sex chromosome must be the same as in the bgen file.**
