# ExtractTransmission

This command line extracts transmitted alleles from hard calls in trio data.


### Requirements

This software requires **java 8** or newer. A 64 bits system is recommended. 


### Installation

The latest release is the folder named `triogen-X.Y.Z` in the [bin folder](https://github.com/mvaudel/trioGen/tree/master/bin), where _X_, _Y_, and _Z_ represent the version number. You can get it by cloning this repository or downloading the folder content. The executable is named accordingly `triogen-X.Y.Z.jar`.


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
-g/--geno          The vcf file.
-f/--fam           The trio identifiers file.
-o/--out           The stem to use for the output files.
```


### Additional Parameters

```
-z/--timeOut       The number of days before timeout, default is 365..
-t/--test          If present, runs only othe first 1000 variants.
```

### Command line example

This example runs from the folder of the cloned repo.

```
java -Xmx16G -cp bin/triogen-0.0.1/triogen-0.0.1.jar no.uib.triogen.cmd.transmission.ExtractTransmission -g /vcfFolder/vcf/22.vcf.gz -f /triadFolder/fulltriads -o /outputFolder/22
```

An example to run a test on the Hunt clound can be found [here](https://github.com/mvaudel/trioGen/blob/master/scripts/test.sh).

### Format requirements

Input files can be gzipped or not. Output files are gzipped. 

- The vcf file should containt the 8 mandatory columns at the beginning. The hard calls of the genotypes of the samples should be indicated in the beginning of the field by a `0` or `1` indicating _ref_ and _alt_ alleles, respectively, and separated by a pipe (`|`) or a slash (`/`). Example: `0|1:bla:bla:bla`. Columns shoud be tab-separated.

- The trio identifiers file should contain three columns with headers containing the identifiers of the children, father, and mother, in that order. Columns should be space-separated.

- For a given stem `/my/stem`, four output files will be created, `/my/stem_h1.gz`, `/my/stem_h2.gz`, `/my/stem_h3.gz`, `/my/stem_h4.gz`, containing _h1_, _h2_, _h3_, and _h4_, respectively. The files contain the eigth mandatory vcf columns followed by the _h_ of each sample. Columns are tab-separated. The columns are ordered alphabetically, the lines are not ordered. **There is no guarantee that the line order corresponds to the vcf file, and no guarantee that the order is consistent between _h1_, _h2_, _h3_, and _h4_.**

> Note: only minimal sanity check is conducted.


### Paths

It is recommended to use absolute paths. The output folder must exist.


### Memory Settings

The processing might require more than the default memory provided to the Java virtual machine. If you encounter memoyr issues please increase the maximum memory setting using the `Xmx` option. 

Example, for a maximum of 64GB of memory:
```
java -Xmx64G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.transmission.ExtractTransmission [parameters]
```


### Errors, questions, and bug report

Errors are reported as stacktrace. Please use the [Issue Tracker]() for questions and bug reports. Please include the stacktrace with any information helping the debugging.



