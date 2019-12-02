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

timix <- function(rho, omega, dat, y, f, r, grp, ftest = F) {

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

  # F-test
  if(ftest == T) {
    n <- nrow(dat)
    p <- ncol(X)
    V <- diag(n) * res$Ve
    for(i in 1:length(Z_list)){
      V <- V + Z_list[[i]] %*% K_list[[i]] %*% t(Z_list[[i]]) *
        res$Vu * res$weights[i]
    }
    sigma_sq <- (t(dat[, y] - X %*% res$betahat) %*% solve(V) %*%
                   (dat[, y] - X %*% res$betahat)) / (n - p)
    p_vec <- testFixed(X, y, res$betahat, sigma_sq, V, n, p)

    # return results
    res <- c(res, p_vec = list(p_vec))
  }
  return(res)
}

#' obtain negative log likelihood of timix function for optimization process
#' @param prm prm[1] is rho and prm[2] is omega.
#' This notification is for oprimization process
#' @export

timixRev <- function(prm, dat, y, f, r, grp) {
  res <- timix(rho = prm[1], omega = prm[2], dat = dat,
               y = y, f = f, r = r, grp = grp)
  return(-1 * res$loglik)
}

#' Mixed model for the synchronized time sereis by optimization.
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

timixCov <- function(dat, y, f, r, grp, grid = T, trace = F){
  if(grid == T){
    prm <- gridprm(rho_vc = seq(0.1, 0.9, 0.2), omega_vc = seq(0.1, 0.9, 0.2),
                   dat = dat, y = y, f = f, r = r, grp = grp)
  }else{
    prm <- c(0, 0)
  }
  if(trace == T){
    print(paste("rho=", prm[1], "omega =", prm[2]))
  }
  opt <- nlminb(start = prm, timixRev,
                dat = dat, y = y, f = f, r = r, grp = grp,
                lower = c(0, 0), upper = c(1 - 1e-5, 1 - 1e-5),
                control = list(trace = trace))
  rho = opt$par[1]; omega = opt$par[2]
  res <- timix(rho, omega, dat, y, f, r, grp, ftest = T)
  return(list(res, opt))
}

#' Mixed model for the synchronized time sereis by grid-search and optimization.
#'
#' @description
#' timixGrd function is a mixed model function to solve synchronized time sereis.
#' The effects of synchronization and serial correlation are treated as the random effects.
#' The variance-covariance matrix of random effects are prescribed by the parameters rho and omega.
#' timixGrd function estimate these two parameters through grid-search and oprimization.
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
#' res <- timixGrd(dat, y = 1, f = 2:4, r = 2:4, grp = 5)
#'

timixGrd <- function(dat, y, f, r, grp, trace = 0){
  prm <- gridprm(rho_vc = seq(0.1, 0.9, 0.2), omega_vc = seq(0.1, 0.9, 0.2),
                 dat = dat, y = y, f = f, r = r, grp = grp)
  if(trace != 0){print(paste("rho=", prm[1], "omega =", prm[2]))}
  opt <- nlminb(start = prm, timixRev,
                dat = dat, y = y, f = f, r = r, grp = grp,
                lower = c(0, 0), upper = c(1 - 1e-6, 1 - 1e-6),
                control = list(trace = trace))
  rho = opt$par[1]; omega = opt$par[2]
  res <- timix(rho, omega, dat, y, f, r, grp, ftest = T)
  return(list(res, opt))
}

#' Grid search of rho and omega
#' @description
#' This function conduct the grid search for selecting the best parameter set of rho and omega
#' @param rho_ve vector of the parameter rho
#' @param omega_ve vector of the parameter omega

gridprm <- function(rho_vc, omega_vc, dat, y, f, r, grp) {
  prm_mt <- expand.grid(rho_vc, omega_vc)
  lik_vc <- numeric(nrow(prm_mt))
  g_elm <- 1
  for(i in 1:nrow(prm_mt)) {
    lik_vc[g_elm] <- timixRev(prm = c(prm_mt[g_elm, 1], prm_mt[g_elm, 2]),
                              dat, y, f, r, grp)
    g_elm <- g_elm + 1
  }
  prm <- c(prm_mt[which.min(lik_vc), 1], prm_mt[which.min(lik_vc), 2])
  return(prm)
}
