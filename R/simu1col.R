#' Function to generate data of one column
#' @param n_time Length of the time series
#' @param n_grp Number of groups
#' @param rho Parameter of serial correlation
#' @param omega Parameter of synchronization
#' @param rndVar Set random numbers for the each group mean.
#' Default is FALSE
#' @param mean mean value (hyper parameter) in the distribution of the group mean.
#' Default is 0
#' @param sd sd value (hyper parameter) in the distribution of the group mean.
#' Default is 0
#'
#' @examples
#' n_time <- 100
#' n_grp <- 5
#'
#' y <- simuOneCol(n_time, n_area, rho = 0.5, omega = 0.5,
#'  rndVar = T, mean = 0, sd = 1)
#' x0 <- rep(1, n_time * n_grp) # Intercept
#' x1 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 1)
#' x2 <- simuOneCol(n_time, n_grp, rho = 0.5, omega = 0)
#' area <- makeGrp(n_time, n_grp)
#'
#' dat <- data.frame(y, x0, x1, x2, area)
#'
#' @export

simuOneCol <- function(n_time, n_grp, rho, omega,
                       rndVar = F, mean = 0, sd = 0) {
  # make a vector to input generated data
  res <- numeric(n_time * n_grp)
  Sig <- makeSyncMat(omega, n_grp)

  # generate initial values
  init_val <- makeRandNum(rndVar, mean, sd, Sig)
  res_elm <- elmSet(n_time, n_grp, 1)
  res[res_elm] <- intm_val <- init_val

  # generate intermediate values
  for(t in 2:n_time) {
    eps <- makeRandNum(F, 0, 0, Sig)
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
makeRandNum <- function(rndVar, mean, sd, Sig){
  if(rndVar){
    m_vec <- rnorm(nrow(Sig), mean, sd)
  }else{
    m_vec <- rep(mean, nrow(Sig))
  }
  rand_num <- mvtnorm::rmvnorm(1, mean = m_vec, sigma = Sig)
  return(rand_num)
}
