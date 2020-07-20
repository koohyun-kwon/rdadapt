#' Inverse Modulus
#'
#' Calculates the inverse modulus for the regression function at a point
#' problem. More specifically, this calcultes
#' \eqn{\omega^{-1}(b, \Lambda_{V+}(C),\Lambda_{V+}(C')) }
#'
#' @param b point where the inverse modulus is evaluated at.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}.
#' @param X n by k design matrix.
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
#' @param C_pair C_pair a pair of smoothness parameters \eqn{(C, C')}.
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
#' minb_fun(C_pair, X, mon_ind)

minb_fun <- function(C_pair, X, mon_ind, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  minb <- min(C_pair[1] * Norm(Vplus(X, mon_ind)) +
                C_pair[2] * Norm(Vminus(X, mon_ind)))

  return(minb)
}
