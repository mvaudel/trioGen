# ExtractTransmission

This command line extracts transmitted alleles from hard calls in trio data.

**Standard command line**

```
java -cp bin/triogen-0.0.1/triogen-0.0.1.jar no.uib.triogen.cmd.transmission.ExtractTransmission [parameters]
```

** Standard parameters **

```
-h/--help          Display help text
-v/--version       Display version
```

** Mandatory Parameters **

```
-g/--geno          The vcf file.
-f/--fam           The trio identifiers file.
-o/--out           The stem to use for the output files.
```

** Format requirements ***

Input files can be gzipped or not. Output files are gzipped. 

- The vcf file should containt the 8 mandatory columns at the beginning. The hard calls of the genotypes of the samples should be indicated in the beginning of the field by a `0` or `1` indicating _ref_ and _alt_ alleles, respectively, and separated by a pipe (`|`) or a slash (`/`). Example: `0|1:bla:bla:bla`. Columns shoud be tab-separated.

- The trio identifiers file should contain three columns with headers containing the identifiers of the children, father, and mother, in that order. Columns should be space-separated.

- For a given stem `/my/stem`, four output files will be created, `/my/stem_h1.gz`, `/my/stem_h2.gz`, `/my/stem_h3.gz`, `/my/stem_h4.gz`, containing _h1_, _h2_, _h3_, and _h4_, respectively. The files contain the eigth mandatory vcf columns followed by the _h_ of each sample. Columns are tab-separated. The columns are ordered alphabetically, but the lines are not ordered. There is no guarantee that the line order corresponds to the vcf file, and no guarantee that the order is consistent between _h1_, _h2_, _h3_, and _h4_.

_Note: only very minimal sanity check is conducted._


** Paths ***

We recommend using absolute paths.
