## MendelianCheck

This command writes a report on the prevalence of Mendelian errors in the given data set.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).

The tested allele is selected as detailed [here](../AlleleFrequency.md).

Mendelian error prevalence is estimated as detailed [here](../MendelianErrors.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.mendelian_check.MendelianCheck [parameters]
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
-f/--fam                  The trio identifiers file.
-o/--out                  The file where to write the results.
```


#### Additional Parameters

```
-vi/--variantId           File listing the variants to include in the analysis. Default: process all variants in the genotypes file.
-af/--afThreshold         Allele frequency threshold. 0.005 excludes all alleles of variants with frequency < 0.5% or > 99.5%. Default: 0.005.
-nv/--nVariants           The number of variants to process in parallel. Default is the number of cores on the machine.
-z/--timeOut              The number of days before timeout, default is 365.
```

### Processing

For each variant, the number of tested alleles for each haplotype is counted and compared to the expected number of alleles based on the allele frequency. 


### Output

The output is a text file with one line per variant.



