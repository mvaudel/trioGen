# ExtractTransmission

This command line extracts transmitted alleles from hard calls in trio data.


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.transmission.ExtractTransmission [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


### Standard parameters

```
-h/--help          Display help text
-v/--version       Display version
```


### Mandatory Parameters

```
-g/--geno          The genotypes file.
-gf/--genoFormat   The genotypes file format. 0: VCF, 1: Sanger VCF.
-gf/--genoFormat   The genotypes file format. 0: VCF, 1: Sanger VCF.
-f/--fam           The trio identifiers file.
-o/--out           The stem to use for the output files.
```


### Additional Parameters

```
-nv/--nVariants    The number of variants to process in parallel. Default is 8.
-z/--timeOut       The number of days before timeout, default is 365.
-t/--test          If present, runs only othe first 1000 variants.
```

### Command line example

This example runs from the folder of the cloned repository.

```
java -Xmx16G -cp bin/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.transmission.ExtractTransmission -g src/main/resources/transmission/test_transmission.vcf -gf 1 -f src/main/resources/transmission/test_trio -o src/test/resources/transmission/result
```


### Output

For a given stem `/my/stem`, four output files will be created, `/my/stem_h1.gz`, `/my/stem_h2.gz`, `/my/stem_h3.gz`, `/my/stem_h4.gz`, containing _h1_, _h2_, _h3_, and _h4_ matrices, respectively. For each line, the variant ID is followed by the _h_ of each sample. Columns are tab-separated. The columns are ordered alphabetically, the lines are not ordered. **There is no guarantee that the line order corresponds to the vcf file, and no guarantee that the order is consistent between _h1_, _h2_, _h3_, and _h4_.**



