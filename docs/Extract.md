# Extract

This command line splits the lines of the linear model output based on the levels of a category column, and/or extracts specific columns.


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.results.Extract [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


### Standard parameters

```
-h/--help          Display help text
-v/--version       Display version
```


### Mandatory Parameters

```
-i/--input         The result file of the LinearModel command. Can be gzipped.
-o/--out           Stem of the file where to write the output. Output will be gzipped.
```


### Additional Parameters

```
-cat/--category      The categories columns as comma separated list, one file will be produced per level. Example: phenotype,variantId. Default: no category.
-val/--value         The columns to include in the results file as comma separated list. Example: variantId,h_B1,h_B1_se,h_B1_p. Default: all columns.
```

### Command line example

The example below runs extracts the `variantId`, `h_B1_p`, `h_B2_p`, `h_B3_p`, and `h_B4_p` columns and makes one output file per phenotype. This can typically be used to extract the values necessary to make Manhattan plots.

```
java -Xmx16G -cp bin/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.results.Extract -i myResults.gz -o myResultsTrimmed -c phenotype -v variantId,h_B1_p,h_B2_p,h_B3_p,h_B4_p
```

