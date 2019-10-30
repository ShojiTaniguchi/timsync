#' Function to generate dependent variable
#'
#' @description
#' genDepnd function is one of the simulation functions of the package timsync.
#' genDepnd generates the value of dependent variable from the independent x and
#' its/their coefficient beta.
#' @param x Independent variable. Vector or matrix.
#' @param beta Coefficients of x.
#' @param rho Parameter of serial correlation.
#' @param omega Parameter of synchronization.
#' @export

genDepnd <- function(x, beta, n_grp, rho, omega) {
  if(class(x) == "matrix") {
    n_time <- nrow(x) / n_grp
    y <- numeric(nrow(x))
  }else if(class(x) == "numeric"){
    n_time <- length(x) / n_grp
    y <- numeric(length(x))
  }
  Sig <- makeSyncMat(omega, n_grp)

  # Generate initial value
  res_elm <- elmSet(n_time, n_grp, 1)
  y_old <- genDepndIni(x, res_elm, beta, Sig, omega)
  y[res_elm] <- y_old

  # Generate values sequencially
  for(i in 2:n_time) {
    res_elm <- elmSet(n_time, n_grp, i)
    y_old <- genDepndSeq(x, y_old, res_elm, beta, Sig, rho, omega)
    y[res_elm] <- y_old
  }
  return(y)
}

#' Function to generate Initial value
#' @param rel_elm Position of the focul element in the data x and y.

genDepndIni <- function(x, res_elm, beta, Sig, omega) {
  if(class(x) == "matrix") {
    y_1 <- x[res_elm] %*% beta
  }
  else if(class(x) == "numeric"){
    y_1 <- x[res_elm] * beta
  }
  else{
    warning("x should be matrix or numeric")
    break
  }
  y_1 <- y_1 + makeRandNum(mean = 0, Sig)
  return(y_1)
}

#' Function to generate dependent variable sequentially
#' @param x independent variables
#' @param y_old previous values one lag before the generating time point

genDepndSeq <- function(x, y_old, res_elm, beta, Sig, rho, omega) {
  if(class(x) == "matrix") {
    y_t <- x[res_elm] %*% beta + y_old * rho + makeRandNum(mean = 0, Sig)
  }else{    #if(class(x) == "numeric")
    y_t <- x[res_elm] * beta + y_old * rho + makeRandNum(mean = 0, Sig)
  }
  return(y_t)
}
