
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/210819589.svg)](https://zenodo.org/badge/latestdoi/210819589)


# TrioGen

TrioGen is a lightweight open-source bioinformatic tool to conduct genome-wide association studies in {child, mother, father} trios. TrioGen is written in java and freely available under the permissive [GNU General Public License v3.0](https://github.com/mvaudel/trioGen/blob/master/LICENSE).


## Getting Started

### Requirements

TrioGen requires **java 8** or newer installed. A 64 bits system is recommended. TrioGen should work on any operating system but Linux is recommended. Documentation, testing, and support are all designed for Linux.


### Installation

The latest release is the folder named `triogen-X.Y.Z` in the [bin folder](https://github.com/mvaudel/trioGen/tree/master/bin), where _X_, _Y_, and _Z_ represent the version number. You can get it by cloning this repository or downloading the folder content. The executable is named accordingly `triogen-X.Y.Z.jar` and should be accompanied by folders names `lib` and `resources`.


## Functionalities

Documentation on the functionalities and associated command lines is available [here](docs/CommandLines.md).


## Files and formats

Documentation on the requires file formats is available [here](docs/FileFormats.md).


## Code

The code is available at the [GitHub repository](https://github.com/mvaudel/trioGen), JavaDoc can be found [here](docs/javadoc/index.html). 


## Performance

General considerations to manage the computational performance of the tool can be found [here](docs/Performance.md).


## Citation

TrioGen does not have an accompanying scientific publication yet. To cite the software please use the URL of the repository and the [DOI](https://zenodo.org/badge/latestdoi/210819589). To cite the association models using the transmitted alleles please cite [Chen _et al._](https://doi.org/10.1101/737106).


## Errors, questions, and bug report

### Code of Conduct

As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Convenant](https://www.contributor-covenant.org/) [Code of Conduct for Open Source Projects](docs/CODE_OF_CONDUCT.md).

Despite our efforts at enforcing good practices in our work, like every software, TrioGen will crash at one time or another, will fail to cover specific use cases, will not perform as well as expected in specific cases. We apologize for any inconvenience and will try to fix things to the best of our capabilities. However, we would like to remind the user that we are scientists and not professional programmers, that we maintain code and provide support on our free time. 

We welcome bug reports, suggestions of improvements, and contributions. Please do not hesitate to [open an issue](https://github.com/mvaudel/trioGen/issues) or a [pull request](https://github.com/mvaudel/trioGen/pulls).


### Bug report

When the tool crashes, an error message with a stacktrace should appear with guidance on how to fix the problem.

Example:
```
java.lang.IllegalArgumentException: The phenotypes file does not contain a "child_id" column.
        at no.uib.triogen.model.pheno.PhenotypesHandler.<init>(PhenotypesHandler.java:114)
        at no.uib.triogen.processing.association.linear_model.LinearModelComputer.run(LinearModelComputer.java:109)
        at no.uib.triogen.cmd.association.LinearModel.run(LinearModel.java:92)
        at no.uib.triogen.cmd.association.LinearModel.main(LinearModel.java:61)

```

If the error messages do not allow you to solve the problem, please use the [Issue Tracker](https://github.com/mvaudel/trioGen/issues) for questions and bug reports. Please include the stacktrace with the information necessary to reproduce the problem.


### Troubleshooting 

General considerations on troubleshooting and commonly encountered errors can be found [here](docs/Troubleshooting.md).

