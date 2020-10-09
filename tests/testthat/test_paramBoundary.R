library(timsync)
load("../cpue_kuroshio_dat.Rdata")
timsync_res <- timixCov(cpue_kuroshio_dat,
                        y = 1, f = c(2, 3, 4), r = c(2, 3, 4), grp = 5, trace = F, grid = F)

mmll <- mmLogLik(rho = timsync_res[[2]]$par[1], omega = timsync_res[[2]]$par[2],
                 dat = cpue_kuroshio_dat, y = 1, f = 2:4, r = 2:4, grp = 5,
                 b = timsync_res[[1]]$betahat,
                 w = timsync_res[[1]]$weights, Ve = timsync_res[[1]]$Ve, Vu = timsync_res[[1]]$Vu, method = "ML")

test_that("Fixed effects", {
  expect_equal(timsync_res[[1]]$betahat[1], -7.226817e-01, tolerance = 0.0001)
  expect_equal(timsync_res[[1]]$betahat[2], 2.084502e-02, tolerance = 0.0001)
  expect_equal(timsync_res[[1]]$betahat[3], -1.363433e-05, tolerance = 0.0001)
})

test_that("Random effects", {
  expect_equal(sqrt(timsync_res[[1]]$Vu * timsync_res[[1]]$weights)[1], 0, tolerance = 1e-5)
  expect_equal(sqrt(timsync_res[[1]]$Vu * timsync_res[[1]]$weights)[2], 0.00751, tolerance = 1e-5)
  expect_equal(sqrt(timsync_res[[1]]$Vu * timsync_res[[1]]$weights)[3], 0.00381, tolerance = 1e-5)
})

test_that("Spatiotemporal correlations", {
  expect_equal(timsync_res[[2]]$par[1], 0.146, tolerance = 1e-3)
  expect_equal(timsync_res[[2]]$par[2], 0.324, tolerance = 1e-3)
})

test_that("LogLikelihood", {
  expect_equal(mmll, -135.4506, tolerance = 1e-4)
})
