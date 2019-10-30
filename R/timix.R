#' Mixed model for the synchronized time sereis with given covariance parameters.
#'
#' @description
#' timix function is a mix model function to solve synchronized time sereis.
#' The effects of synchronization and serial correlation are treated as the random effects.
#' The variance-covariance matrix of random effects are prescribed by the parameters rho and omega.
#'
#' @param rho parameter for serial correlation
#' @param omega parameter for synchronization
#' @param dat data for the analysis
#' @param y column number of the dependent variable
#' @param f column number of the fixed effects
#' @param r column number of the random effects
#' @param grp column number of the group
#'
#' @examples
#' n_time <- 50
#' n_grp <- 5

#' y <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 0.5, mean = 0)
#' x0 <- rep(1, n_time * n_grp) # Intercept
#' x1 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 1, mean = 0)
#' x2 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 0, mean = 0)
#' area <- makeGrp(n_time, n_grp)

#' dat <- data.frame(y, x0, x1, x2, area)
#' timix(rho = 0.5, omega = 0.5, dat, y = 1, f = 2:4, r = 2:4, grp = 5)
#' @export

timix <- function(rho, omega, dat, y, f, r, grp) {

  # number of groups and time points
  n_group <- nlevels(dat[, grp])
  n_time <- nrow(dat) / n_group

  # R matrix for serial correlation
  R <- matrix(1, nr = n_time, nc = n_time)
  for ( i in 1:nrow(R) ) {
    for ( j in 1:ncol(R) ) {
      R[i, j] <- rho ^ abs(i-j)
    }
  }

  # Omega_0 matrix for synchronization
  Omega <- matrix(1, nr = n_group, nc = n_group)
  for ( i in 1:nrow(Omega) ) {
    for ( j in 1:ncol(Omega) ) {
      if (i != j) { Omega[i, j] <- omega }
    }
  }
  Omeaga_0 <- Omega %x% R

  # make design matrix: X
  X <- dat[, f]
  X <- as.matrix(X)
  colnames(X) <- colnames(dat)[f]

  # make design matrix: Z
  Z_list <- NULL
  for(i in r){
    Z_list <- c(Z_list, list(diag(dat[, i])))
  }

  # make variance-covariance matrix
  K_list <- NULL
  for(i in r){
    K_list <- c(K_list, list(Omeaga_0))
  }

  # fit model
  res <- tsEMMREMLMultiKernel(y = dat[, y], X = X, Zlist = Z_list, Klist = K_list,
                                     varbetahat = T, varuhat = T, test = T)
  return(res)
}

#' obtain negative log likelihood of timix function for optimization process
#' @param prm prm[1] is rho and prm[2] is omega.
#' This notification is for oprimization process

timixRev <- function(prm, dat, y, f, r, grp) {
  res <- timix(rho = prm[1], omega = prm[2], dat = dat,
               y = y, f = f, r = r, grp = grp)
  return(-1 * res$loglik)
}

#' Mixed model for the synchronized time sereis by estimating covariance parameters.
#'
#' @description
#' timixCov function is a mixed model function to solve synchronized time sereis.
#' The effects of synchronization and serial correlation are treated as the random effects.
#' The variance-covariance matrix of random effects are prescribed by the parameters rho and omega.
#' timixCov function estimate these two parameters through oprimization.
#'
#' @param dat data for the analysis
#' @param y column number of the dependent variable
#' @param f column number of the fixed effects
#' @param r column number of the random effects
#' @param grp column number of the group
#' @examples
#' n_time <- 50
#' n_grp <- 5
#' y <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 0.5, mean = 0)
#' x0 <- rep(1, n_time * n_grp) # Intercept
#' x1 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 1, mean = 0)
#' x2 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 0, mean = 0)
#' area <- makeGrp(n_time, n_grp)
#' dat <- data.frame(y, x0, x1, x2, area)
#' res <- timixCov(dat, y = 1, f = 2:4, r = 2:4, grp = 5)
#'
#' @export

timixCov <- function(dat, y, f, r, grp){
  opt <- try(nlminb(start = c(0, 0), timixRev,
                    dat = dat, y = y, f = f, r = r, grp = grp,
                    lower = c(0, 0), upper = c(1, 1)))
  rho = opt$par[1]; omega = opt$par[2]
  res <- timix(rho, omega, dat, y, f, r, grp)
  return(list(res, opt))
}
