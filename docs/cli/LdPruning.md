## LdPruning

This command LD-prunes result files.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.ld_pruning.LdPruning [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-res/results              The results file to prune. Can be any text file, gzipped or not.
-l/--ldMatrix             The ld matrix file as generated using the [_LdMatrix_](LdMatrix.md) command. If LD matrix files are computed per chromosome, replace the chromosome name with '{contig}'.
-f/--fam                  The trio identifiers file.
-o/--out                  The file where to write the results.
```


#### Additional Parameters

```
-r/--minR2                The minimal ld r2 to consider two markers in LD. Default: 0.05. Min value: value used when generating the LD matrix file.
-p/--maxP                 The maximal p-value to consider. Default: 1e-6.
-id/--idColName           The name of the variant id column. Default: 'variantId'.
-pn/--pColName            The name of the p-value column. Default: 'h.intercept.p'.
-cn/--contigColName       The name of the contig column. Default: 'contig'.
-phn/--phenoColName       The name of the phenotype column. Ignored if not provided.
-s/--separator            Separator for the columns. Default: '\t'.
```

### Processing

For each phenotype, variants are selected from lowest p-values until reaching `maxP` excluding variants with an LD R2 >= `minR2` with any of the already selected variants.


### Output

The entire lines of the variants passing the pruning are written to the output as they were in the input.



