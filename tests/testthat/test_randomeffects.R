library(timsync)
library(nlme)

set.seed(100)
x0 <- rep(1, 100)
x1 <- rnorm(100, 0, 1)
Omega <- makeOmega(n_time = 20, n_grp = 5, rho = 1, omega = 0)
y <- genDat(x = x1, beta = 1, n_time = 20, n_grp = 5, z = x1,
            OmegaMat = list(Omega), sig = c(1, 1), rho = 0, omega = 0)
grp <- rep(paste(1:5, "g", sep = ""), each = 20)
dat <- data.frame(y = y$objVar, x0, x1, grp)
rslt_lme <- lme(y ~ x1, dat, random = ~ x1 - 1 | grp, method = "REML")
rslt_lme_out <- capture.output(rslt_lme)
pos <- grep("StdDev", rslt_lme_out)
stdDev <- rslt_lme_out[pos]
stdDev <- unlist(strsplit(stdDev, split = " "))[-1]
stdDev1 <- as.numeric(stdDev)
rslt_sum <- summary(rslt_lme)
b_nlme1 <- rslt_sum$tTable[, 1]

timsync_res <- timix(rho = 1, omega = 0, dat, y = 1, f = 2:3, r = 3, grp = 4)
Vu1 <- timsync_res$Vu
Ve1 <- timsync_res$Ve
b_timsync1 <- timsync_res$betahat

set.seed(10)
x0 <- rep(1, 100)
x1 <- rnorm(100, 0, 1)
Omega <- makeOmega(n_time = 20, n_grp = 5, rho = 1, omega = 0)
y <- genDat(x = x1, beta = 1, n_time = 20, n_grp = 5, z = x1,
            OmegaMat = list(Omega), sig = c(1, 1), rho = 0, omega = 0)
grp <- rep(paste(1:5, "g", sep = ""), each = 20)
dat <- data.frame(y = y$objVar, x0, x1, grp)
rslt_lme <- lme(y ~ x1, dat, random = ~ x1 - 1 | grp, method = "REML")
rslt_lme_out <- capture.output(rslt_lme)
pos <- grep("StdDev", rslt_lme_out)
stdDev <- rslt_lme_out[pos]
stdDev <- unlist(strsplit(stdDev, split = " "))[-1]
stdDev2 <- as.numeric(stdDev)
rslt_sum <- summary(rslt_lme)
b_nlme2 <- rslt_sum$tTable[, 1]

timsync_res <- timix(rho = 1, omega = 0, dat, y = 1, f = 2:3, r = 3, grp = 4)
Vu2 <- timsync_res$Vu
Ve2 <- timsync_res$Ve
b_timsync2 <- timsync_res$betahat

test_that("Variance are the same between timsync and nlme", {
  expect_equal(Vu1, stdDev1[1]^2, tolerance = 10e-3)
  expect_equal(Ve1, stdDev1[2]^2, tolerance = 10e-3)
  expect_equal(Vu2, stdDev2[1]^2, tolerance = 10e-3)
  expect_equal(Ve2, stdDev2[2]^2, tolerance = 10e-3)
})

test_that("Coefficients are the same between timsync and nlme", {
  expect_equal(as.numeric(b_timsync1), as.numeric(b_nlme1), tolerance = 10e-3)
  expect_equal(as.numeric(b_timsync2), as.numeric(b_nlme2), tolerance = 10e-3)
})

