#' F-test for the fixed effects
#' @param X design matrix for the fixed effects.
#' @param y dependent variable.
#' @param beta a vector of the estimated coefficients.
#' @param sigma_sq the estimated sigma square, a component of the denominator of the F value.
#' @param V the variance-covariance matrix of the dependent variavble.
#' @param n the number of data
#' @param p the number of fixed effects

testFixed <- function(X, y, beta, sigma_sq, V, n, p) {
  p_vec <- numeric(p)
  for(i in 1:p){
    H <- matrix(0, nrow = p, ncol = 1)
    H[i, 1] <- 1
    fval <- t(t(H) %*% beta) %*%
      solve(t(H) %*% solve(t(X) %*% solve(V) %*% X) %*% H) %*%
      t(H) %*% beta / sigma_sq
    p_vec[i] <- pf(fval, 1, n - p, lower.tail = F)
  }
  return(p_vec)
}
