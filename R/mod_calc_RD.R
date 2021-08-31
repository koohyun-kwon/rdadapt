#' Inverse Modulus
#'
#' Calculates the inverse modulus for the regression function at a point
#' problem. More specifically, this calcultes
#' \eqn{\omega^{-1}(b, \Lambda_{V+}(C),\Lambda_{V+}(C')) }.
#'
#' @param b point where the inverse modulus is evaluated at.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}.
#' @param X n by d design matrix.
#' @param mon_ind index of the monotone variables.
#' @param sigma standard deviation of the error term (either length 1 or n).
#' @param swap indicator for whether we take (C', C) instead of (C, C').
#'
#' @return the value of inverse modulus given \eqn{\omega = b}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' C_pair <- c(0.5, 1)
#' sigma <- c(1, 2, 1)
#' invmod(1, C_pair, X, mon_ind, sigma)
#' invmod(1, c(Inf, 1), X, mon_ind, sigma)
invmod <- function(b, C_pair, X, mon_ind, sigma, swap = FALSE){

  if (!(length(sigma) %in% c(1, nrow(X)))) {
    stop("sigma must have length 1 or n")
  }

  K <- b * K_fun(b, C_pair, X, mon_ind, swap)
  res <- sqrt(sum(K^2 / sigma^2))

  return(res)
}

#' Minimum Modulus
#'
#' Calculates the smallest possible modulus value, i.e., \eqn{\omega(0)}.
#'
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}.
#' @param X A data matrix.
#' @param mon_ind index number for monotone variables.
#' @param swap indicator for whether we take (C', C) instead of (C, C').
#'
#' @return the value of the smallest possible modulus value.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' C_pair <- c(0.5, 1)
#' C_pair2 <- c(Inf, 1)
#' minb_fun(C_pair, X, mon_ind)
#' minb_fun(C_pair2, X, mon_ind)
minb_fun <- function(C_pair, X, mon_ind, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  b1 <- C_pair[1] * Norm(Vplus(X, mon_ind))
  b2 <- C_pair[2] * Norm(Vminus(X, mon_ind))

  # Adjustment for the case where C = Inf
  b1[Norm(Vplus(X, mon_ind)) == 0] = 0
  b2[Norm(Vminus(X, mon_ind)) == 0] = 0

  minb <- min(b1 + b2)

  return(minb)
}

#' Inverse Modulus for RD Parameter
#'
#' Calculates the inverse modulus for the RDD problem.
#'
#' @param b point where the inverse modulus is evaluated at.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}.
#' @param Xt \eqn{n_t} by \eqn{d} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{d} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#' @param swap indicator for whether we take (C', C) instead of (C, C').
#'
#' @return A list of four components. \code{bt} gives the modulus value
#' corresponding to the treated observations, while \code{delta_t} gives
#' the value of the inverse modulus corresponding to the treated observations.
#' \code{bc} and \code{delta_c} give the analogous values for the control observations.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' C_pair <- c(0.1, 0.2)
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' invmod_RD(1, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' C_pair[1] = Inf
#' invmod_RD(1, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)
invmod_RD <- function(b, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  ## Derivative of the square of the minization problem
  deriv_bt <- function(bt) {

    bc <- b - bt

    om_inv_t_der <- bt * K_fun(bt, C_pair, Xt, mon_ind) / sigma_t^2
    om_inv_t_der <- sum(om_inv_t_der)

    om_inv_c_der <- bc * K_fun(bc, C_pair, Xc, mon_ind, swap = TRUE) / sigma_c^2
    om_inv_c_der <- sum(om_inv_c_der)

    return(om_inv_t_der - om_inv_c_der)
  }

  minbt <- minb_fun(C_pair, Xt, mon_ind)
  minbc <- minb_fun(C_pair, Xc, mon_ind, swap = T)
  minb <- minbt + minbc

  if(b == minb){  # delta = 0 case

    bt <- minbt

  }else if(deriv_bt(minbt) * deriv_bt(b - minbc) > 0){ # Corner solution

    bt <- b - minbc

  }else{

    bt_sol <- stats::uniroot(deriv_bt, c(minbt, b - minbc), tol = .Machine$double.eps^10)
    bt <- bt_sol$root
  }

  delta_t <- invmod(bt, C_pair, Xt, mon_ind, sigma_t)
  delta_c <- invmod(b - bt, C_pair, Xc, mon_ind, sigma_c, swap = TRUE)

  res <- list(bt = bt, delta_t = delta_t, bc = b - bt, delta_c = delta_c)

  return(res)

}

#' Modulus of Continuity for the RD Parameter
#'
#' Calculates the modulus of continuity for the RD parameter.
#'
#' @param delta a nonnegative value of \eqn{\delta}.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}.
#' @param Xt \eqn{n_t} by \eqn{k} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{k} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#' @param swap indicator for whether we take (C', C) instead of (C, C').
#'
#' @return A list of four components. \code{bt} gives the modulus value
#' corresponding to the treated observations, while \code{delta_t} gives
#' the value of the inverse modulus corresponding to the treated observations.
#' \code{bc} and \code{delta_c} give the analogous values for the control observations.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' C_pair <- c(0.1, 0.2)
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' modsol_RD(1, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c, swap = FALSE)
#' C_pair[1] <- Inf
#' modsol_RD(1, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c, swap = FALSE)
modsol_RD <- function(delta, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  minbt <- minb_fun(C_pair, Xt, mon_ind)
  minbc <- minb_fun(C_pair, Xc, mon_ind, swap = TRUE)
  minb <- minbt + minbc

  eqn_fun <- function(b) {

    invmod_RD_res <- invmod_RD(b, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)
    delta_t <- invmod_RD_res$delta_t
    delta_c <- invmod_RD_res$delta_c

    res <- sqrt(delta_t^2 + delta_c^2) - delta
    return(res)
  }

  maxint <- 100 # An arbitrary large number; doesn't affect the result
  solve <- stats::uniroot(eqn_fun, c(minb, maxint), extendInt = "upX",
                          tol = .Machine$double.eps^10)
  bsol <- solve$root

  res <- invmod_RD(bsol, C_pair, Xt, Xc, mon_ind, sigma_t, sigma_c)

  return(res)
}
