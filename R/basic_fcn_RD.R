#' Positive Part
#'
#' Calculates component-wise \code{max(X,0)} for a matrix X.
#'
#' @param X A data matrix.
#'
#' @return A matrix in the same dimension of \code{X}, after \code{max(x,0)}
#' was applied to each component \code{x} of \code{X}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' pos(X)
pos <- function(X) {

  return(pmax(X, 0))
}

#' Function for \eqn{(x)_{V+}}
#'
#' Calculates rowwise \eqn{(x)_{V+}} for each row \eqn{x}
#' of a matrix \code{X}.
#'
#' @param X A data matrix.
#' @param mon_ind index number for monotone variables.
#'
#' @return A matrix in the same dimension of \code{X},
#' after \eqn{(x)_{V+}} was applied to each row \eqn{x} of \code{X}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' Vplus(X, mon_ind)
Vplus <- function(X, mon_ind) {

  Xplus <- X
  Xplus[, mon_ind] <- pos(X[, mon_ind])

  return(Xplus)
}

#' Function for \eqn{(x)_{V-}}
#'
#' Calculates rowwise \eqn{(x)_{V-}} for each row \eqn{x}
#' of a matrix \code{X}.
#'
#' @param X A data matrix.
#' @param mon_ind index number for monotone variables.
#'
#' @return A matrix in the same dimension of \code{X},
#' after \eqn{(x)_{V-}} was applied to each row \eqn{x} of \code{X}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' Vminus(X, mon_ind)
Vminus <- function(X, mon_ind) {

  return(Vplus(-X, mon_ind))
}

#' (Weighted) \eqn{L_p} Norm
#'
#' Calculates rowwise \eqn{||x||_p} for each row \eqn{x}
#' of a matrix \code{X}.
#'
#' \code{X} should be a matrix, not a vector.
#'
#' @param X A numeric matrix.
#' @param p Order of the norm.
#' @param invw Inverse weights for each component,
#' with the length equal to either 1 or \code{ncol(X)}.
#'
#' @return A vector of the dimension equal to the \code{nrow(X)}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' Norm(X)
#' weights <- c(1, 2)
#' Norm(X = X, invw = weights)
Norm <- function(X, p = 2, invw = 1) {

  n <- nrow(X)
  d <- ncol(X)

  invw <- matrix(rep(invw, each = n), nrow = n, ncol = d)
  Xw <- X / invw

  return(rowSums(abs(Xw)^p)^(1/p))
}

#' Kernel Function for Adaptive Procedures
#'
#' Calculates \eqn{K(x/b; C, C')}.
#'
#' @param b "bandwidth" parameter.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')}
#' @param X data matrix.
#' @param mon_ind index number for monotone variables.
#' @param swap we take (C', C) instead of (C, C') if \code{swap = TRUE}.
#'
#' @return A numeric vector with the length equal to nrow(X);
#' if b = 0, returns 0 vector.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' C_pair <- c(0.5, 1)
#' K_fun(1, C_pair, X, mon_ind)
K_fun <- function(b, C_pair, X, mon_ind, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  if(b == 0){

    res <- rep(0, nrow(X))

  }else{

    res <- pos(1 - (C_pair[1] / b) * Norm(Vplus(X, mon_ind)) -
               (C_pair[2] / b) * Norm(Vminus(X, mon_ind)))
  }

  return(res)
}
