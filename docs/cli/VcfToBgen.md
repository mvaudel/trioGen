## VcfToBgen

This command converts a phased vcf file to bgen v. 1.3 with phased haplotypes.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.vcf_to_bgen.VcfToBgen [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-i/--input                 The vcf file.
-o/--output                The bgen file.
```


### Processing

Since genotyping probabilities are not phased, only hard-called genotypes are exported. Pobabilities are coded using an 8 bits precision. Genotypes are compressed using Zstd.


### Output

A bgen file and an index file containing the location of every variant.



