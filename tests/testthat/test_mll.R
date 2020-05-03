library(timsync)

set.seed(100)
x0 <- rep(1, 100)
x1 <- rnorm(100, 0, 1)
Omega <- makeOmega(n_time = 20, n_grp = 5, rho = 1, omega = 0)
y <- genDat(x = x1, beta = 1, n_time = 20, n_grp = 5, z = x1,
            OmegaMat = list(Omega), sig = c(1, 1), rho = 0, omega = 0)
grp <- rep(paste(1:5, "g", sep = ""), each = 20)
dat <- data.frame(y = y$objVar, x0, x1, grp)
library(nlme)
rslt_lme <- lme(y ~ x1, dat, random = ~ x1 - 1 | grp, method = "ML")
rslt_lme_out <- capture.output(rslt_lme)
pos <- grep("StdDev", rslt_lme_out)
stdDev <- rslt_lme_out[pos]
stdDev <- unlist(strsplit(stdDev, split = " "))[-1]
stdDev <- as.numeric(stdDev)
rslt_sum <- summary(rslt_lme)
nlme_mll <- rslt_sum$logLik

timsync_mll <- mmLogLik(rho = 1, omega = 0, dat = dat, y = 1, f = 2:3, r = 3, grp = 4,
                        b = rslt_sum$tTable[, 1],
                        w = 1, Vu = stdDev[1]^2, Ve = stdDev[2]^2,
                        method = "ML")

test_that("Maximum likelihood are the same between timsync and nlme", {
  expect_equal(nlme_mll, timsync_mll, tolerance = 10e-3)
})
