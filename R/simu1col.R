#' Function to generate data of one column
#' @param n_time Length of the time series
#' @param n_grp Number of groups
#' @param rho Parameter of serial correlation
#' @param omega Parameter of synchronization
#'
#' @examples
#' n_time <- 100
#' n_grp <- 5
#'
#' y <- simuOneCol(n_time, n_area, rho = 0.5, omega = 0.5, mean = 0)
#' x0 <- rep(1, n_time * n_area) # Intercept
#' x1 <- simuOneCol(n_time, n_area, rho = 0.5, omega = 1, mean = 0)
#' x2 <- simuOneCol(n_time, n_area, rho = 0.5, omega = 0, mean = 0)
#' area <- makeGrp(n_time, n_area)
#'
#' dat <- data.frame(y, x0, x1, x2, area)
#'
#' @export

simuOneCol <- function(n_time, n_grp, rho, omega, mean = 0) {
  # make a vector to input generated data
  res <- numeric(n_time * n_grp)
  Sig <- makeSyncMat(omega, n_grp)

  # generate initial values
  init_val <- makeRandNum(mean, Sig)
  res_elm <- elmSet(n_time, n_grp, 1)
  res[res_elm] <- intm_val <- init_val

  # generate intermediate values
  for(t in 2:n_time) {
    eps <- makeRandNum(mean, Sig)
    intm_val <- intm_val * rho + eps
    res_elm <- elmSet(n_time, n_grp, t)
    res[res_elm] <- intm_val
  }
  return(res)
}

#' Function to set element numbers of the intermediate values at time i
#'@param n_time Length of the time series
#'@param n_grp Number of areas
#'@param i Time point i
elmSet <- function(n_time, n_grp, t) {
  (0:(n_grp - 1)) * n_time + t
}

#' Function to  make group
#'@param n_time Length of the time series
#'@param n_grp Number of areas
#'@export
makeGrp <- function(n_time, n_grp) {
  rep(paste(1:n_grp, "area", sep = ""), each = n_time)
}

#' Make the variance-covariance matrix to prescribe synchronization
#' @param omega Parameter of synchronization
#' @param n_grp Number of areas
#' @export
makeSyncMat <- function(omega, n_grp) {
  Sig <- matrix(omega, nrow = n_grp, ncol = n_grp)
  diag(Sig) <- 1
  return(Sig)
}

#' Make random numbers with synchronization
#' @param mean Mean of the random numbers
#' @param Sig Variance-Covariance matrix to prescribe synchronization
makeRandNum <- function(mean, Sig){
  rand_num <- mvtnorm::rmvnorm(1, mean = rep(mean, nrow(Sig)), sigma = Sig)
  return(rand_num)
}
