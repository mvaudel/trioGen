
##
#
# This script generates dummy data to test covariates handling.
#
##


covariatesDF <- data.frame(
    cov1 = c(1, 1, 0, 0, 1),
    cov2 = c(1, 0, 1, 0, 0),
    cov3 = c(0, 1, 0, 1, 0)
)

svdResults <- svd(
    x = covariatesDF,
    nu = 5
)
svdResults

x <- c(1, 2, 1, 2, 3)
y <- c(2, 3, 2, 3, 4)

x0 <- c(1, 2, 3)
y0 <- c(2, 3, 4)

lm(y0 ~ x0)

t(svdResults$u) %*% x
t(svdResults$u) %*% y


uDF <- read.table("tmp/svd/z_bmi0_u.gz", header = F, stringsAsFactors = F)
uMatrix <- as.matrix(uDF)
tu <- t(uMatrix)
nrow(tu)
ncol(tu)

sDF <- read.table("tmp/svd/z_bmi0_s.gz", header = F, stringsAsFactors = F)
sMatrix <- as.matrix(sDF)

vDF <- read.table("tmp/svd/z_bmi0_v.gz", header = F, stringsAsFactors = F)
vMatrix <- as.matrix(vDF)


rawPheno <- t(as.matrix(read.table("tmp/svd/z_bmi0_rawPheno.gz", header = F, stringsAsFactors = F)))
product <- t(as.matrix(read.table("tmp/svd/z_bmi0_covariateBaseValues.gz", header = F, stringsAsFactors = F)))
result <- t(as.matrix(read.table("tmp/svd/z_bmi0_result.gz", header = F, stringsAsFactors = F)))

