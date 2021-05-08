## LocusZoom

This command extracts the data needed to make locus zoom plots from the results of the [_LinearModel_](LinearModel.md) command in text tables.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).

A connection to the Ensembl API must be available for gene mapping.


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.locus.LocusZoom [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-res/results              The results file of the [_LinearModel_](LinearModel.md) command.
-l/--ldMatrix             The LD matrix file as generated using the [_LdMatrix_](LdMatrix.md) command. If LD matrix files are computed per chromosome, replace the chromosome name with '{contig}'.
-f/--fam                  The trio identifiers file.
-vi/--variantId           File listing the variants to query.
-o/--out                  The file where to write the data needed to build a locus zoom plot.
```


#### Additional Parameters

```
-pn/--phenoName           The name of the phenotype to export the locus zoom data on. Example: pheno1. Default: all phenotypes found for the given variants.
-d/--dist                 The maximum distance in bp to consider around a variant. Should be below the distance used to create the ld matrix and, if gene mapping is conducted, below 2.5e6. Default: 1000000.
-g/--geneCoordinates      The file where to write the gene coordinates. If none provided gene mapping will be skipped.
-b/--buildNumber          The number of the build to use if gene mapping is conducted. e.g. 38 for GRCh38. Default: 38.
-log/--log                The file where to write the log. Default: next to the output.
```

### Processing

For each phenotype, variants are selected from lowest p-values until reaching `maxP` excluding variants with an LD R2 >= `minR2` with any of the already selected variants.


### Output

The entire lines of the variants passing the pruning are written to the output as they were in the input.



