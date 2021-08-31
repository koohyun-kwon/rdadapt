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

  Xplus <- X # initialization
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
#' @param p Order of the norm; the default is \code{p = 1} (L1 norm).
#' @param invw Inverse weights for each component,
#' with the length equal to either 1 or \code{ncol(X)}; the default is
#' \code{invw = 1}.
#'
#' @return A vector of the dimension equal to the \code{nrow(X)}.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' Norm(X)
#' weights <- c(1, 2)
#' Norm(X = X, invw = weights)
Norm <- function(X, p = 1, invw = 1) {

  n <- nrow(X)
  d <- ncol(X)

  # Error message for wrong weights
  if(length(invw) != 1){
    if(length(invw) != d){
      stop("dimension of invw should be either 1 or ncol(X)")
    }
  }

  invw <- matrix(rep(invw, each = n), nrow = n, ncol = d)
  Xw <- X / invw

  return((rowSums(abs(Xw)^p))^(1/p))
}


#' Kernel Function for Adaptive Procedures
#'
#' Calculates \eqn{K(x/b; C, C')} in our previous draft.
#'
#' @param b a scalar "modulus of continuity" parameter;
#' can take the value of 0.
#' @param C_pair a pair of smoothness parameters \eqn{(C, C')};
#' can take the values of \code{Inf}.
#' @param X data matrix.
#' @param mon_ind index number for monotone variables, with length
#' \code{ncol(X)}.
#' @param swap we take (C', C) instead of (C, C') if \code{swap = TRUE}.
#'
#' @return A numeric vector with the length equal to \code{nrow(X)};
#' if b = 0, returns 0 vector.
#' @export
#'
#' @examples X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
#' mon_ind <- c(1, 2)
#' C_pair <- c(0.5, 1)
#' K_fun(1, C_pair, X, mon_ind)
#' K_fun(1, c(Inf, 1), X, mon_ind)
K_fun <- function(b, C_pair, X, mon_ind, swap = FALSE){

  if (swap) {
    C_pair <- C_pair[2:1]
  }

  if(b == 0){
    res <- rep(0, nrow(X))
  }else{
    comp1 <- (C_pair[1] / b) * Norm(Vplus(X, mon_ind))
    comp2 <- (C_pair[2] / b) * Norm(Vminus(X, mon_ind))

    # Adjustment for C = Inf
    comp1[Norm(Vplus(X, mon_ind)) == 0] = 0
    comp2[Norm(Vminus(X, mon_ind)) == 0] = 0

    res <- pos(1 - comp1 - comp2)
  }

  return(res)
}


#' Row-wise Minimum
#'
#' Calculates row-wise minimum of a matrix by looping over columns.
#'
#' This function is useful when \code{nrow} is much greater than \code{ncol}.
#'
#' @param mat_var input matrix
#'
#' @return a vecvor with a dimension \code{nrow(mat_var)}.
#' @export
#'
pmin_mat <- function(mat_var){

  J <- ncol(mat_var)

  if(J == 1){

    res <- mat_var

  }else{

    res <- mat_var[, 1]

    for(j in 2:J){

      res <- pmin(res, mat_var[, j])
    }
  }
  return(res)
}


#' Standard Deviation Given Weights
#'
#' Calculate a standard deviation of a linear estimator given weights and conditional
#' standard deviations.
#'
#'
#' @param w_t weights for treated observations.
#' @param w_c weights for control observations.
#' @param sigma_t conditional standard deviations for the treated.
#' @param sigma_c conditional standard deviations for the control.
#'
#' @return a scalar value.
#' @export
#'
#' @examples w_t <- rep(1/10, 10)
#' w_c <- rep(1/5, 5)
#' sigma_t <- rep(1/2, 10)
#' sigma_c <- rep(1/3, 5)
#' sd_w(w_t, w_c, sigma_t, sigma_c)
sd_w <- function(w_t, w_c, sigma_t, sigma_c){

  res2 <- sum(w_t^2 * sigma_t^2) + sum(w_c^2 * sigma_c^2)
  res <- sqrt(res2)

  return(res)
}
