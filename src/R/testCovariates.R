
##
#
# This script generates dummy data to test covariates handling.
#
##


covariatesDF <- data.frame(
    cov1 = c(1, 1, 0, 0, 0),
    cov2 = c(1, 0, 0, 0, 0),
    cov3 = c(0, 1, 0, 0, 0)
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
