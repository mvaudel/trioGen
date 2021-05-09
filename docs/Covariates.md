## Covariates

Covariates should be numeric. For factors, we recommend [one-hot endcoding](https://en.wikipedia.org/wiki/One-hot). These covariates are typically PCs or batches. We only conduct linear standardization, non-linear operations must be conducted when normalizing the phenotypes, see [gamlss](https://www.gamlss.com/) for examples. Please note that this will introduce a difference in how genotypes and phenotypes are standardized.

For a given phenotype, if a covariate is missing, it will be imputed using the median of the 200 covariates with nearest phenotypic value.

Similarly as in [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/), the covariates are then projected out for both the genotypes and phenotypes using the [SingularValueDecomposition implementation](https://commons.apache.org/proper/commons-math/javadocs/api-3.6/org/apache/commons/math3/linear/SingularValueDecomposition.html) of the [Commons Math library](http://commons.apache.org/proper/commons-math/).

