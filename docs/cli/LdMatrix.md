## LdMatrix

This command computes Ld matrices from bgen files.


### General considerations

Files need to be formatted as detailed [here](../FileFormats.md).


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.ld_matrix.LdMatrix [parameters]
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
-d/--dist                 The maximum distance in bp to consider around a variant. Default: 500000.
-r/--minR2                The minimal ld r2 to report (inclusive). Default: 1e-6.
-af/--afThreshold         Allele frequency threshold. 0.001 excludes all alleles of variants with frequency < 0.1% or > 99.9%. Only variants with at least two alleles passing the threshold will be inspected. Default: 0.001.
-nv/--nVariants           The number of variants to process in parallel. Default is 8.
-z/--timeOut              The number of days before timeout, default is 365.
```

### Processing

All alleles of all variants within the given bp window are compared. LD R2 is computed using the parents in the trio file, when no genotype is found for the parents of a trio, the genotype of the child is used.


### Performance considerations

Depending on the window size and the number of samples, RAM requirement can become important possibly limiting the number of chromosomes that can be run in parallel. For ~100,000 samples and a max distance of 500,000 kb, this command requires approximately 20 GB of RAM. If the number of CPU used is lower than the number of variants to process in parallel, this is due to the reading and parsing of the bgen file, consider using ssd discs and removing low maf variants if this happens. 

Please make sure that your disk is large enough to save the results. Depending on the threshold, the file size is around 1-10 kB per variant.


### Output

The output file LD values between all the alleles of all the variants considered in a binary block-compressed format as detailed [here](../FileFormats.md).



