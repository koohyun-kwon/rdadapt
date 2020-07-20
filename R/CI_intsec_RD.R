#' Covariance Between Two Estimators
#'
#' Calculates the covariance between two estimators,
#' \eqn{Lhat_j(\delta)} and \eqn{Lhat_k(\delta)}.
#'
#' @param delta a nonegative scalar value;
#' it can be left unspecified if
#' (\code{ht_j}, \code{ht_k}, \code{hc_j}, \code{hc_k}) are specified.
#' @param Cj the smoothness parameter corresponding to the first estimator.
#' @param Ck the smoothness parameter corresponding to the second estimator.
#' @param Cbar the largest smoothness parameter.
#' @param Xt \eqn{n_t} by \eqn{k} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{k} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#' @param ht_j the modulus value for the treated observations,
#' corresponding to the first smoothness parameter;
#' it can be left unspecified if \code{delta} is specified.
#' @param hc_j the modulus value for the contorl observations,
#' corresponding to the first smoothness parameter;
#' it can be left unspecified if \code{delta} is specified.
#' @param ht_k the modulus value for the treated observations,
#' corresponding to the second smoothness parameter;
#' it can be left unspecified if \code{delta} is specified.
#' @param hc_k the modulus value for the control observations,
#' corresponding to the second smoothness parameter;
#' it can be left unspecified if \code{delta} is specified.
#'
#' @return a scalar covariance value.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] > 0 & X[, 2] > 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' cov_calc(1, 0.2, 0.4, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
cov_calc <- function(delta, Cj, Ck, Cbar, Xt, Xc, mon_ind,
                     sigma_t, sigma_c, ht_j, hc_j, ht_k, hc_k){

  if(missing(ht_j) | missing(hc_j) | missing(ht_k) | missing(hc_k)){

    hres_j <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_j <- hres_j$ht
    hc_j <- hres_j$hc

    hres_k <- bw_adpt(delta, Ck, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_k <- hres_k$ht
    hc_k <- hres_k$hc
  }

  num_tj <- K_fun(b = ht_j, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind)
  num_tk <- K_fun(b = ht_k, C_pair = c(Cbar, Ck), X = Xt, mon_ind = mon_ind)
  num_t <- sum(num_tj * num_tk / sigma_t^2)
  sum_t <- ht_j * ht_k * num_t / delta^2

  num_cj <- K_fun(b = hc_j, C_pair = c(Cj, Cbar), X = Xc, mon_ind = mon_ind)
  num_ck <- K_fun(b = hc_k, C_pair = c(Ck, Cbar), X = Xc, mon_ind = mon_ind)
  num_c <- sum(num_cj * num_ck / sigma_c^2)
  sum_c <- hc_j * hc_k * num_c / delta^2

  res <- sum_t + sum_c
  return(res)
}

#' Covariance Matrix Calculation for Multiple Estimators
#'
#' Calculates the covariance matrix of \eqn{J} estimators,
#' \eqn{Lhat_1(\delta)},..., \eqn{Lhat_J(\delta)}.
#'
#' @param delta a nonegative scalar value;
#' it can be left unspecified if
#' \code{hmat} is specified.
#' @param Cvec a sequence of smoothness parameters
#' @inheritParams cov_calc
#' @param hmat a \eqn{J} by {2} matrix of modulus values;
#' it can be left unspecified if \code{delta} is specified.
#'
#' @return a \eqn{J} by \eqn{J} covariance matrix.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] > 0 & X[, 2] > 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' cov_mat_calc(1, (1:5)/5, Xt, Xc, mon_ind, sigma_t, sigma_c)
cov_mat_calc <- function(delta, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, hmat){

  J <- length(Cvec)
  Cbar <- max(Cvec)

  if(missing(hmat)){

    hmat <- matrix(0, nrow = J, ncol = 2)

    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$ht
      hmat[j, 2] <- hres_j$hc
    }
  }

  res <- diag(1, J, J)

  for(j in 1:(J - 1)){

    for(k in (j + 1):J){

      Cj <- Cvec[j]
      Ck <- Cvec[k]
      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]
      ht_k <- hmat[k, 1]
      hc_k <- hmat[k, 2]

      res[j, k] <- cov_calc(delta, Cj, Ck, Cbar, Xt, Xc, mon_ind,
                            sigma_t, sigma_c, ht_j, hc_j, ht_k, hc_k)
      res[k, j] <- res[j, k]
    }
  }

  return(res)
}


#' Quantile of Maximum of Multivariate Normal Random Variable
#'
#' Calculates the quantile of maximum of a multivariate normal random variable,
#' by a Monte Carlo simulation method.
#'
#' @param covmat covariance matrix of the multivariate normal random variable.
#' @param alpha desired upper quantile value.
#' @param num_sim number of simulations used to calculate the quantile;
#' the default is \code{10^4}.
#'
#' @return \eqn{(1 - \alpha)}th quantile of maximum of a multivariate normal random variable.
#' @export
#'
#' @examples covmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' max_Q(covmat, 0.05)
max_Q <- function(covmat, alpha, num_sim = 10^4){

  J <- nrow(covmat)
  rn <- MASS::mvrnorm(n = num_sim, mu = rep(0, J), Sigma = covmat)
  rn_max <- apply(rn, 1, max)

  res <- stats::quantile(rn_max, 1 - alpha)
  return(res)
}
