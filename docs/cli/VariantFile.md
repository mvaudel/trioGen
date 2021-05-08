## VariantFile

This command makes a variant file from a list of identifiers.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).

A connection to the Ensembl API must be available.


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.variant_file.VariantFile [parameters]
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
-vi/--variantId           A file containing a list of variant rsids to query, one line per variant.
-s/--source               A value to write in the column 'source' for future reference.
-o/--out                  The file where to write the results.
```


#### Additional Parameters

```
-b/--build                The build to use when querying Ensembl as a number 37: grch37, 38: grch38. Default: 37.
-p/--population           The reference population to use as available in Ensembl. See https://rest.ensembl.org/documentation/info/variation_populations for details. Default: 1000GENOMES:phase_3:GBR.
-r/--r2                   The minimal r2 to allow for a proxy. Default: 0.2.
```

### Processing

If a variant id is missing, a proxy is queried in Ensembl. The ids of variants where no proxy is found is available in the log.


### Output

A variant file with one line per variant.



