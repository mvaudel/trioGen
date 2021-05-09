## VariantData

This command extracts the standardized genotypes and phenotypes for a set of variants in a text table. This allows convenient follow-up analyses like interaction or MR studies.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).

Covariates are handled as detailed [here](../Covariates.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.variant_data.VariantData [parameters]
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
-vi/--variantId           File listing the variants to include in the analysis along with the weights.
-p/--phenoFile            The phenotypes and covariates file.
-pn/--phenoName           List of the names of the phenotypes in the phenotype file (Example: pheno1,pheno2).
-f/--fam                  The trio identifiers file.
-o/--out                  The file where to write the results.
```


#### Additional Parameters

```
-cg/--covariate_general   List of the names of the covariates to use for all phenotypes. Example: pc1,pc2,pc3,pc4,pc5,pc6.
-cs/--covariate_specific  File containing the names of the covariates to use for a specific phenotype in json format.
-id/--childId             The name of the column containing the child id. Default: child_SentrixID.
-vl/--variantLog          If present, writes a log for every variant next to the results file.
```

### Processing

For each variant, genotypes and phenotypes are adjusted for covariates and the matrix that is otherwise used for the [_LinearModel_](LinearModel.md) command is exported to a text file. 


### Output

The output is a text file with one line per variant and trio, with the haplotypes and phenotypes in column. Note that sample identifiers are not included.



