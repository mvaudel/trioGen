## LdValue

This command returns variants in LD with a given variant.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.ld_value.LdValue [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-l/--ldMatrix             The LD matrix file as generated using the [_LdMatrix_](LdMatrix.md) command. If LD matrix files are computed per chromosome, replace the chromosome name with '{contig}'.
-vi/--variantId           File listing the variants to query.
-o/--out                  The file where to write the results.
```


### Output

A json file containing all variants in the LD matrix file for which an LD value is present. LD R2 is reported for all combinations of alleles available.



