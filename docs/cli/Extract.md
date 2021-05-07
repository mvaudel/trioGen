## Extract

This command line splits the lines of the [_LinearModel_](cli/LinearModels.md) command output based on the levels of a category column, and/or extracts specific lines or columns.


### Command line

```
java -Xmx16G -cp your/folder/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.extract.Extract [parameters]
```

> Note: you need to replace `your/folder` by the folder where the release is installed, and `Z.Y.Z` by the version number.


#### Standard parameters

```
-h/--help                 Display help text
-v/--version              Display version
```


#### Mandatory Parameters

```
-i/--input                The result file of the [_LinearModel_](cli/LinearModels.md) command.
-o/--out                  Stem of the file where to write the output. Output will be gzipped and indexed.
```


#### Additional Parameters

```
-sv/--split_by_variant    Splits the output creating one file per variant. Default: no split.
-sp/--split_by_pheno      Splits the output creating one file per phenotype. Default: no split.
-col/--columns            The columns containing association summary statistics to include in the results file as comma separated list. Example: cmf.Bc,cmf.Bc.se,cmf.Bc.p. Default: all columns.
-id/--variantId           The ids of the variants to include in the results file as comma separated list. Example: rs123,rs456,rs789. Default: all variants.
-p/--pheno                The phenotypes to include in the results file as comma separated list. Example: pheno1,pheno2,pheno3. Default: all phenotypes.
```

### Command line example

The example below runs extracts the `variantId`, `cmf.Bc.p`, `cmf.Bm.p`, `cmf.Bf.p`, `cmf.h.p`, and `cmf.intercept.p` columns and makes one output file per phenotype. This can typically be used to extract the values necessary to make Manhattan plots.

```
java -Xmx16G -cp bin/triogen-X.Y.Z/triogen-X.Y.Z.jar no.uib.triogen.cmd.extract.Extract -i myResults.gz -o myResultsTrimmed -c phenotype -v variantId,cmf.Bc.p,cmf.Bm.p,cmf.Bf.p,cmf.h.p,cmf.intercept.p
```

