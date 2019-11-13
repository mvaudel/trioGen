# TrioGen

TrioGen is a lightweight open-source bioinformatic tool to conduct genetic analyses in trios. TrioGen is written in java and freely available under the permissive [GNU General Public License v3.0](https://github.com/mvaudel/trioGen/blob/master/LICENSE).

TrioGen uses genotypes from trios defined as {child, father, mother} to identify transmitted and non-transmitted alleles and run association with phenotypes. Throughout the code and documentation, the nomenclature introduced by [Chen _et al._](https://doi.org/10.1101/737106) is used as illustrated below. 

![h](docs/illustrations/h.png?raw=true "Nomenclature by Chen et al.")

For a given variant, TrioGen uses phased genotypes to extract the alleles transmitted by the mother to the child (Am1) or non-transmitted (Am2), and the alleles transmitted by the father to the child (Af1) or non-transmitted (Af2). The number of alternative alleles transmitted or not by each parent is noted _h_ in the model by [Chen _et al._](https://doi.org/10.1101/737106), as detailed in the table below.

| h | Description |
| ---- | ----------- |
| h1 | Number of alternative allele transmitted to the child from the mother |
| h2 | Number of alternative allele of the mother not transmitted to the child |
| h3 | Number of alternative allele transmitted to the child from the father |
| h4 | Number of alternative allele of the father not transmitted to the child |

## Command Lines

1. ExtractTransmission: This command line extracts _h_ matrices. [[Documentation]](docs/ExtractTransmission.md)

2. LineaModel: This command line runs simple linear regression on phenotypes against _h_. [[Documentation]](docs/ExtractTransmission.md)


## Getting Started

### Requirements

TrioGen requires **java 8** or newer installed. A 64 bits system is recommended.


### Installation

The latest release is the folder named `triogen-X.Y.Z` in the [bin folder](https://github.com/mvaudel/trioGen/tree/master/bin), where _X_, _Y_, and _Z_ represent the version number. You can get it by cloning this repository or downloading the folder content. The executable is named accordingly `triogen-X.Y.Z.jar`.


## Files and Formats

> **Warning: only minimal sanity check is conducted.**

### Paths

It is recommended to use absolute paths. The output folder must exist.


### Input and Output

Input files can be gzipped or not. Output files are gzipped. 

> **Warning: there is not guarantee that the order of the lines or columns in the output is the same as in the input.**

### Genotyping files

Currently, only monoallelic variants from diploid chromosomes are supported. Each variant must have a unique identifier. Genotypes of the children must be phased, with the `A|B` notation representing the allele inherited from the mother to the left (`A`) and from the father to the right (`B`). Genotyping files must be in a valid [_vcf_ format](https://www.htslib.org/doc/vcf.html). All samples should be in the same vcf files. One vcf file is processed at a time, if genotypes are spread among multiple vcf files, _e.g._ one per chromosome, use multiple command lines.

Two parsers are available and can be selected when running the command lines (see details on each command line).

1. VCF: This option uses the generic [htslib](htslib.org) parser and can parse [all vcf files supported by htslib](https://github.com/samtools/hts-specs). The vcf file must be [compressed and indexed](https://genometoolbox.blogspot.com/2014/09/how-to-index-vcf-file.html).

2. Sanger VCF: This option uses a custom parser tested on VCF files produced by the Sanger Imputatoin service only, allowing a faster processing of the file but is not compatible with all vcf files. This parser is provided with no guarantee, and we recommend using the _Generic VCF_ parser for testing and if this parser does not work on your files. The vcf file must containt the 8 mandatory columns at the beginning. The hard calls of the genotypes of the samples must be indicated in the beginning of the field by a `0` or `1` indicating _ref_ and _alt_ alleles, respectively, and separated by a pipe (`|`) or a slash (`/`). Example: `0|1:bla:bla:bla`. Columns shoud be tab-separated. The vcf file can be gzipped or not.


### Trio file

Trios must be provided as text file. The file can be gzipped or not. The file must contain three columns with headers containing the identifiers of the children, father, and mother, in that order. Columns must be space-separated.


## Performance

### I/O

According to our preliminary tests, I/O is the main limitation in terms of speed. You can improve this by using SSD discs.

### Threading

You can set the number of variants to process in parallel using the command line options. Note, however, that if the performance on your setup is limited by the I/O, this will not have much of an effect. Note also that some variant-specific tasks are parallelized when possible using all available resources. You can override the number of threads used by default using the `-Djava.util.concurrent.ForkJoinPool.common.parallelism` argument. 

Example with 32 threads:
```
java -Djava.util.concurrent.ForkJoinPool.common.parallelism=32 -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.Command [parameters]
```

It is possible to run multiple instances of TrioGen from the same _jar_ file. The preformance gain of processing multiple vcf files in parallel depends on your setup but should be limited.


### Memory Settings

The processing might require more than the default memory provided to the Java virtual machine. If you encounter memoyr issues please increase the maximum memory setting using the `Xmx` option. 

Example, for a maximum of 64 GB of memory:
```
java -Xmx64G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.Command [parameters]
```


## Errors, questions, and bug report

Errors are reported as stacktrace. Please use the [Issue Tracker](https://github.com/mvaudel/trioGen/issues) for questions and bug reports. Please include the stacktrace with the information necessary to reproduce the problem.

