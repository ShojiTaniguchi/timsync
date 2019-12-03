#' Mixed model for the synchronized time sereis by optimization, when omega = 0
#'
#' @description
#' timixCovRho function is a mixed model function to solve synchronized time sereis when omega = 0.
#' The effects of serial correlation is treated as the random effects.
#' The variance-covariance matrix of random effects are prescribed by the parameter rho.
#' timixCovRho function estimate the parameter through oprimization.
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
#' res <- timixCovRho(dat, y = 1, f = 2:4, r = 2:4, grp = 5)
#'
#' @export

timixCovRho <- function(dat, y, f, r, grp, grid = T, trace = F){
  if(grid == T){
    prm <- gridprm(rho_vc = seq(0.1, 0.9, 0.2), omega_vc = 0,
                   dat = dat, y = y, f = f, r = r, grp = grp)
  }else{
    prm <- c(0, 0)
  }
  if(trace == T){
    print(paste("rho=", prm[1], "omega =", prm[2]))
  }
  opt <- nlminb(start = prm[1], timixRevRho,
                dat = dat, y = y, f = f, r = r, grp = grp,
                lower = 0, upper = 1 - 1e-5,
                control = list(trace = trace))
  rho = opt$par[1]
  res <- timix(rho, omega = 0, dat, y, f, r, grp, ftest = T)
  return(list(res, opt))
}

#' obtain negative log likelihood of timix function for optimization process
#' @param prm prm[1] is rho, while omega is set to be 0.
#' This notification is for oprimization process
#' @export

timixRevRho <- function(prm, dat, y, f, r, grp) {
  res <- timix(rho = prm, omega = 0, dat = dat,
               y = y, f = f, r = r, grp = grp)
  return(-1 * res$loglik)
}
