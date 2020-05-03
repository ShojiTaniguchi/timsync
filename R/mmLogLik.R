#' Function to calculate maximum likelihood
#' @description
#' mmLogLik function is to calculate the maximum likelihood values
#' under the estimated parameters.
#' The appropriate usage is to calculate the ML under the REML parameters
#' estiamted by the timixCov or timix function.

#' @param rho parameter for serial correlation
#' @param omega parameter for synchronization
#' @param dat data for the analysis
#' @param y column number of the dependent variable
#' @param f column number of the fixed effects
#' @param r column number of the random effects
#' @param grp column number of the group
#' @param w weight vector for the random effects variance
#' @param Ve Variance of the residuals
#' @param Vu Variance of the random effects
#' @param method "ML" or "REML"

#' @references
#' Hamazaki, K., & Iwata, H. (2020).
#' RAINBOW: Haplotype-based genome-wide association study
#'  using a novel SNP-set method. PLoS computational biology, 16(2), e1007663.

# mmLogLik <- function(y, X, b, Zlist, Klist, w, Ve, Vu, method = "ML"){
mmLogLik <- function(rho, omega, dat, y, f, r, grp, w, b, Ve, Vu, method){

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
  Omega_0 <- Omega %x% R

  # make y
  y <- dat[, y]

  # make design matrix: X
  X <- dat[, f]
  X <- as.matrix(X)
  colnames(X) <- colnames(dat)[f]
  n <- nrow(X)

  # make design matrix: Z
  Z_list <- NULL
  for(i in r){
    Z_list <- c(Z_list, list(diag(dat[, i])))
  }

  # make variance-covariance matrix
  K_list <- NULL
  for(i in r){
    K_list <- c(K_list, list(Omega_0))
  }

  if(class(X) == "matrix"){
    p <- ncol(X)
  }else if(class(X) == "numeric"){
    p <- 1
  }
  Ks <- 0
  for(i in 1:length(Z_list)){
    Ks <- Ks + Z_list[[i]] %*% K_list[[i]] %*% t(Z_list[[i]]) * w[i]
  }
  H <- Ks + Ve/Vu * diag(n)
  detH <- determinant(H, logarithm = T)
  detH <- as.numeric(detH[[1]])

  LogLik <- 1/2 * (-n * log(2 * pi * Vu) - detH -
                     1/Vu * t(y - X %*% b) %*% solve(H) %*% (y - X %*% b))
  if(method == "REML"){
    LogLik <- LogLik + 1/2 * (p * log(2 * pi * Vu) + log(det(t(X) %*% X)) - log(det(t(X) %*% solve(H) %*% X)))
  }
  return(as.numeric(LogLik))
}
