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
#' max_Q2(covmat, 0.05)
max_Q2 <- function(covmat, alpha, num_sim = 10^4){

  J <- nrow(covmat)
  rn <- MASS::mvrnorm(n = num_sim, mu = rep(0, J), Sigma = covmat)
  rn_max <- apply(rn, 1, max)

  res <- as.numeric(stats::quantile(rn_max, 1 - alpha))
  return(res)
}


#' Quantile of Maximum of Multivariate Normal Random Variable 2
#'
#' Calculates the quantile of maximum of a multivariate normal random variable,
#' using multivariate normal distribution function calculation.
#'
#' @param covmat covariance matrix of the multivariate normal random variable.
#' @param alpha desired upper quantile value.
#'
#' @return \eqn{(1 - \alpha)}th quantile of maximum of a multivariate normal random variable.
#' @export
#'
#' @examples covmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' max_Q(covmat, 0.05)
max_Q <- function(covmat, alpha){

  qres <- mvtnorm::qmvnorm(p = 1 - alpha, tail = "lower.tail", sigma = covmat,
                           algorithm = mvtnorm::Miwa())

  res <- qres$quantile
  return(res)
}

#' Coverage Probability Calculation for Intersection CI
#'
#' Calculates the (non-)coverage probabilities of individual confidence intervals
#' that ensure the proper coverage of the intersection confidence interval.
#'
#' @inheritParams cov_mat_calc
#' @inheritParams max_Q
#' @param delta_init the value of \eqn{\delta} to be used in simulating the quantile;
#' theoretically, its value does not matter asymptotically. Its default value is 1.96.
#' @param hmat_init the matrix of modulus values corresponding to \code{delta_init}
#' and \code{Cvec}; it can be left unspecified.
#'
#' @return a list of two values, \code{del_sol}, the proper value of \eqn{\delta} to be used
#' for individual CIs, and \code{tau_sol}, the non-coverage probability associated with the
#' value of \eqn{\delta}.
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
#' tau_calc((1:5)/5, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
tau_calc <- function(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                    delta_init = 1.96, hmat_init){

  J <- length(Cvec)
  Cbar <- max(Cvec)

  if(missing(hmat_init)){

    hmat_init <- matrix(0, nrow = J, ncol = 2)

    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(delta_init, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat_init[j, 1] <- hres_j$ht
      hmat_init[j, 2] <- hres_j$hc
    }
  }

  covmat <- cov_mat_calc(delta_init, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, hmat_init)
  del_sol <- max_Q(covmat, alpha)
  tau_sol <- 1 - stats::pnorm(del_sol)

  res <- list(del_sol = del_sol, tau_sol = tau_sol)
  return(res)
}

#' Lower Adaptive Confidence Interval
#'
#' Calculates the value of the lower end of the adaptive confidence interval.
#'
#' @inheritParams tau_calc
#' @param Yt outcome value for the treated group observations.
#' @param Yc outcome value for the control group observations.
#' @param hmat the matrix of modulus values corresponding to
#' the optimal \eqn{\delta} and \code{Cvec}; it can be left unspecified.
#'
#' @return a list components \code{CI} which gives the value of
#' the lower end of the confidence interval, \code{CI_vec} which gives the values of
#' all CIs used to compute the adaptive CI, and \code{tau_sol} which gives the value of
#' non-coverage probability used to compute individual CIs.
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
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' CI_adpt_L((1:5)/5, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05)
CI_adpt_L <- function(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                    delta_init = 1.96, hmat_init, hmat){

  J <- length(Cvec)
  Cbar <- max(Cvec)

  if(missing(hmat_init)){

    hmat_init <- matrix(0, nrow = J, ncol = 2)

    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(delta_init, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat_init[j, 1] <- hres_j$ht
      hmat_init[j, 2] <- hres_j$hc
    }
  }

  tau_res <- tau_calc(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                      delta_init, hmat_init)
  tau_sol <- tau_res$tau_sol
  del_sol <- tau_res$del_sol

  resvec <- numeric(J)

  for(j in 1:J){

    Cj <- Cvec[j]

    if(missing(hmat)){

      hres_j <- bw_adpt(del_sol, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      ht_j <- hres_j$ht
      hc_j <- hres_j$hc

    }else{

      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]

    }

    resvec[j] <- c_hat_lower_RD(del_sol, Cj, Cbar, Xt, Xc, mon_ind,
                                sigma_t, sigma_c, Yt, Yc, tau_sol,
                                ht_j, hc_j)
  }

  res_list <- list(CI = max(resvec), CI_vec = resvec, tau_sol = tau_sol)
  return(res_list)
}


#' Title
#'
#' @inheritParams CI_adpt_L
#' @param lower calculate a lower one-sided confidence interval if \code{TRUE};
#' calculate a two-sided CI otherwise.
#' @param hmat_init_L the matrix of modulus values corresponding to \code{delta_init}
#' and \code{Cvec}; it can be left unspecified.
#' @param hmat_L the matrix of modulus values corresponding to
#' the optimal \eqn{\delta} and \code{Cvec}; it can be left unspecified.
#'
#' @return confidence interval endpoints.
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
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' CI_adpt((1:5)/5, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = FALSE)
#' CI_adpt((1:5)/5, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = TRUE)
CI_adpt <- function(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, lower = FALSE,
                    delta_init = 1.96, hmat_init_L, hmat_L){

  if(class(Xt) != "matrix" | class(Xc) != "matrix"){
    stop("Xt and Xc should be matrices")
  }

  if(lower == T){

    res_L <- CI_adpt_L(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                     delta_init = 1.96, hmat_init_L, hmat_L)
    res <- c(res_L$CI, Inf)
  }else{

    Cbar <- max(Cvec)

    res_L <- CI_adpt_L(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha/2,
                       delta_init = 1.96, hmat_init_L, hmat_L)
    CI_L <- res_L$CI

    CI_U <- -c_hat_lower_RD(stats::qnorm(1 - alpha/2), Cbar, Cbar, Xc, Xt, mon_ind,
                            sigma_c, sigma_t, Yc, Yt, alpha/2)

    res <- c(CI_L, CI_U)
  }

  return(res)
}
