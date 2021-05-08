## SimpleScore

This command computes a simple score for each trio based on the sum of reference effect sizes.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.simple_score.SimpleScore [parameters]
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
-id/--childId             The name of the column containing the child id. Default: child_SentrixID.
-z/--timeOut              The number of days before timeout, default is 365.
-vl/--variantLog          If present, writes a log for every variant next to the results file.
```

### Processing

For each haplotype, a risk score is computed by summing the weight of the variants provided. If a variant is not found in the genotypes file, it is ignored. 


### Output

The output is a text file with one line per trio and scores for every haplotype.



