#' Function to generate data
#'
#' @description
#' genDat function is the simulation function of the package timsync.
#' genDat generates any variable for the time series data with synchrony.
#'
#' @param x Independent variable. Vector or matrix.
#' @param beta Coefficients of x.
#' @param n_time Length of the time series
#' @param n_grp Number of groups.
#' @param z Random variables. Vector or matrix.
#' @param OmegaMat List of Variance-Covariance Matrix (Omega matrix) for the random effects.
#' @param sig Vector of the standard deviation of the random variables or error term.
#' The last element should be that of the error term.
#' @param rho Parameter of serial correlation of the error term.
#' @param omega Parameter of synchronization of the error term.
#' @examples
#' ## Generate independent variable.
#' y0 <- genDat(x = 0, beta = 0, n_time = 20, n_grp = 3, rho = 0.3, omega = 0.8)
#' grp <- makeGrp(n_time = 20, n_grp = 3)
#' dat <- data.frame(y0 = y0$objVar, time = rep(1:20, 3), grp = grp)
#' ggplot(dat, aes(x = time, y = y0, color = grp)) + geom_line()
#'
#' ## Generate independent variable with random group means.
#' Omega <- makeOmega(n_time = 20, n_grp = 3, rho = 1, omega = 0)
#' y0 <- genDat(x = 0, beta = 0, n_time = 20, n_grp = 3,
#'  z = 1, OmegaMat = list(Omega), sig = c(1, 1),
#'  rho = 0.3, omega = 0.8)
#' grp <- makeGrp(n_time = 20, n_grp = 3)
#' dat <- data.frame(y0 = y0$objVar, time = rep(1:20, 3), grp = grp)
#' ggplot(dat, aes(x = time, y = y0, color = grp)) + geom_line()
#'
#' ## Generate dependent variable with random effects.
#' x1 <- rnorm(20 * 3, 0, 1)
#' Omega_0 <- makeOmega(n_time = 20, n_grp = 3, rho = 0.5, omega = 0.5)
#' Omega_1 <- makeOmega(n_time = 20, n_grp = 3, rho = 0.5, omega = 0.5)
#' Z <- cbind(rep(1, 20 * 3), x1)
#' y0 <- genDat(x = x1, beta = 5, n_time = 20, n_grp = 3,
#'  z = Z, OmegaMat = list(Omega_0, Omega_1), sig = c(1, 1, 1),
#'  rho = 0.3, omega = 0.8)
#'
#' @export

genDat <- function(x, beta, n_time, n_grp,
                   z = NULL, OmegaMat = NULL, sig = 1, rho, omega){

  n_b <- length(OmegaMat)
  n_s <- length(sig)
  if(n_b + 1 != n_s){
    stop("Length of OmegaMat + 1 should be the length of sig")
  }

  # make a vector to input generated data
  res <- numeric(n_time * n_grp)
  Sig <- makeSyncMat(omega, n_grp)

  # make the random effects variables
  b_lis <- NULL
  rnd_val <- 0
  if(n_b > 0){
    for(i in 1:n_b){
      b_vec <- mvtnorm::rmvnorm(1, mean = rep(0, n_time * n_grp),
                                sigma = sig[i] * OmegaMat[[i]])
      if(n_b == 1){
        rnd_val <- rnd_val + z * b_vec
      }else{
        rnd_val <- rnd_val + z[, i] * b_vec
      }
      b_lis <- c(b_lis, list(b_vec))
    }
  }
  # Make the error term
  Omega_e <- makeOmega(n_time, n_grp, rho, omega)
  err_term <- mvtnorm::rmvnorm(1, mean = rep(0, n_time * n_grp),
                               sigma = sig[[n_s]] * Omega_e)

  if(length(beta) == 1){
    obj <- x * beta + as.vector(rnd_val) + as.vector(err_term)
  }else{
    obj <- as.vector(x %*% beta) + as.vector(rnd_val) + as.vector(err_term)
  }
  return(list(objVar = obj, rndEff = b_lis))
}


#' Function to generate Omega matrix
#' @param n_time Length of the time series
#' @param n_grp Number of groups
#' @param rho Parameter of serial correlation
#' @param omega Parameter of synchronization
#' @export

makeOmega <- function(n_time, n_grp, rho, omega) {
  ## R matrix for serial correlation
  R <- matrix(1, nr = n_time, nc = n_time)
  for ( i in 1:nrow(R) ) {
    for ( j in 1:ncol(R) ) {
      R[i, j] <- rho ^ abs(i-j)
    }
  }

  ## O matrix for synchronization
  O <- matrix(1, nr = n_grp, nc = n_grp)
  for ( i in 1:nrow(O) ) {
    for ( j in 1:ncol(O) ) {
      if (i != j) { O[i, j] <- omega }
    }
  }
  return(O %x% R)
}
