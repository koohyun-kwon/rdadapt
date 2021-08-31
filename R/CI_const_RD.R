#' Between Modulus Calculation
#'
#' Calculates the between moduli \eqn{h_{jt}} and \eqn{h_{jc}}
#' for a given value of \eqn{\delta} and \eqn{C_j}.
#'
#' @param delta a nonnegative value of \eqn{\delta}.
#' @param Cj the smoothness parameter aiming to adapt to.
#' @param Cbar the largest smoothness parameter.
#' @param Xt \eqn{n_t} by \eqn{k} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{k} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#'
#' @return a list of two values, \code{h_jt} and \code{h_jc}.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' bw_adpt(1, 1/2, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' bw_adpt(1, 1/2, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
bw_adpt <- function(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c){

  C_pair = c(Cbar, Cj)

  modres <- modsol_RD(delta, C_pair, Xt, Xc, mon_ind,
                      sigma_t, sigma_c, swap = FALSE)

  ht <- modres$bt
  hc <- modres$bc
  res <- list(ht = ht, hc = hc)

  return(res)
}

#' Estimator for the Adaptive CI
#'
#' Calculates \eqn{Lhat_j(\delta)} in the paper
#'
#' @param delta a nonegative scalar value:
#' it can be left unspecified if \code{ht} and \code{hc} are specified.
#' @param Cj the smoothness parameter aiming to adapt to.
#' @param Cbar the largest smoothness parameter.
#' @param Xt \eqn{n_t} by \eqn{k} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{k} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#' @param Yt outcome value for the treated group observations.
#' @param Yc outcome value for the control group observations.
#' @param ht the modulus value for the treated observations;
#' it can be left unspecified if \code{delta} is specified.
#' @param hc the modulus value for the control observations;
#' it can be left unspecified if \code{delta} is specified.
#' @param ret.w returns weights vector if \code{TRUE}.
#'
#' @return a scalar value of the estimator
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' Lhat_fun_RD (1, 1/2, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc)
#' Lhat_fun_RD (1, 1/2, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc)
Lhat_fun_RD <- function(delta, Cj, Cbar, Xt, Xc, mon_ind,
                        sigma_t, sigma_c, Yt, Yc, ht, hc,
                        ret.w = FALSE) {

  if(missing(ht) | missing(hc)){

    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }

  num_it <- K_fun(b = ht, C_pair = c(Cbar, Cj), Xt, mon_ind) * Yt / sigma_t^2
  denom_it <- K_fun(b = ht, C_pair = c(Cbar, Cj), Xt, mon_ind) / sigma_t^2
  res_t <- sum(num_it) / sum(denom_it)

  num_ic <- K_fun(b = hc, C_pair = c(Cj, Cbar), Xc, mon_ind) * Yc / sigma_c^2
  denom_ic <- K_fun(b = hc, C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c^2
  res_c <- sum(num_ic) / sum(denom_ic)

  if(ret.w == TRUE){

    w_t <- (K_fun(b = ht, C_pair = c(Cbar, Cj), Xt, mon_ind) / sigma_t^2) / sum(denom_it)
    w_c <- (K_fun(b = hc, C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c^2) / sum(denom_ic)
    res <- list(est = res_t - res_c, w_t = w_t, w_c = w_c)
  }else{

    res <- res_t - res_c
  }

  return(res)
}


#' Helper Function for Worst-case Bias Calculation
#'
#' Calculates 0.5 * \eqn{(a_{t} - a_{c})} in the worst-case formula.
#'
#' This correponds to the component in the expression (16) of our paper.
#'
#' @param delta a nonegative scalar value:
#' it can be left unspecified if \code{ht} and \code{hc} are specified.
#' @param Cj the smoothness parameter aiming to adapt to.
#' @param Cbar the largest smoothness parameter.
#' @param Xt \eqn{n_t} by \eqn{k} design matrix for the treated units.
#' @param Xc \eqn{n_c} by \eqn{k} design matrix for the control units.
#' @param mon_ind index number for monotone variables.
#' @param sigma_t standard deviation of the error term for the treated units
#' (either length 1 or \eqn{n_t}).
#' @param sigma_c standard deviation of the error term for the control units
#' (either length 1 or \eqn{n_c}).
#' @param ht the modulus value for the treated observations;
#' it can be left unspecified if \code{delta} is specified.
#' @param hc the modulus value for the control observations;
#' it can be left unspecified if \code{delta} is specified.
#'
#' @return a scalar value
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' a_fun(1, 1/2, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' a_fun(1, 1/2, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
a_fun <- function(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, ht, hc){

  if(missing(ht) | missing(hc)){

    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }

  num_it1 <- K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2
  num_it2 <- Cbar * Norm(Vplus(Xt, mon_ind)) - Cj * Norm(Vminus(Xt, mon_ind))
  num_it <- num_it1 * num_it2
  denom_it <- K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2
  res_t <- sum(num_it) / sum(denom_it)

  num_ic1 <- K_fun(b = hc, C_pair = c(Cj, Cbar), X = Xc, mon_ind = mon_ind) / sigma_c^2
  num_ic2 <- Cj * Norm(Vplus(Xc, mon_ind)) - Cbar * Norm(Vminus(Xc, mon_ind))
  num_ic <- num_ic1 * num_ic2
  denom_ic <- K_fun(b = hc, C_pair = c(Cj, Cbar), Xc, mon_ind) / sigma_c^2
  res_c <- sum(num_ic) / sum(denom_ic)

  res <- 0.5 * (res_t - res_c)
  return(res)

}

#' Standard Deviation of Estimator
#'
#' Calculates the standard deviation of the estimator
#' \eqn{Lhat_j(\delta)}.
#'
#' This corresponds to \eqn{s} function defined in Section 2.2 of our paper.
#'
#' @inheritParams a_fun
#'
#' @return a scalar value
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' sd_Lhat_RD(1, 0.5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' sd_Lhat_RD(1, 0.5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
sd_Lhat_RD <- function(delta, Cj, Cbar, Xt, Xc, mon_ind,
                       sigma_t, sigma_c, ht, hc){

  if(missing(ht) | missing(hc)){

    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }

  res <- (delta / ht) /
    sum(K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2)

  return(res)
}

#' Worst-case Bias
#'
#' Calculates the worst case bias of the estimator
#' \eqn{Lhat_j(\delta)}.
#'
#' This corresponds to expression (16) of our paper
#'
#' @inheritParams a_fun
#'
#' @return a scalar value
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' sup_bias_Lhat_RD(1, 1/2, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' sup_bias_Lhat_RD(1, 1/2, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
sup_bias_Lhat_RD <- function(delta, Cj, Cbar, Xt, Xc, mon_ind,
                             sigma_t, sigma_c, ht, hc){

  if(missing(ht) | missing(hc)){

    hres <- bw_adpt(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht <- hres$ht
    hc <- hres$hc
  }

  res1 <- 0.5 * (ht + hc)
  res2 <- a_fun(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, ht, hc)
  # Same as -0.5 * delta * sd_Lhat_RD
  res3 <- -0.5 * (delta^2 / ht) /
    sum(K_fun(b = ht, C_pair = c(Cbar, Cj), X = Xt, mon_ind = mon_ind) / sigma_t^2)

  res <- res1 + res2 + res3
  return(res)
}

#' Lower Adaptive CI for the RD Parameter
#'
#' Calculates the lower end point of the adaptive CI for the RD parameter.
#'
#' This corresponds to the lower CI defined in Section 4.2 of our paper.
#'
#' @inheritParams a_fun
#' @param Yt outcome value for the treated group observations.
#' @param Yc outcome value for the control group observations.
#' @param tau desired level of non-coverage probability.
#'
#' @return the value of the lower end point of the adaptive CI.
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' c_hat_lower_RD(1, 1/2, 1, Xt, Xc, mon_ind, sigma_t, sigma_c,
#' Yt, Yc, 0.05)
#' c_hat_lower_RD(1, 1/2, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c,
#' Yt, Yc, 0.05)
c_hat_lower_RD <- function(delta, Cj, Cbar, Xt, Xc, mon_ind,
                           sigma_t, sigma_c, Yt, Yc, tau, ht, hc) {

  lhat <- Lhat_fun_RD(delta, Cj, Cbar, Xt, Xc, mon_ind,
                      sigma_t, sigma_c, Yt, Yc, ht, hc)

  sup_bias <- sup_bias_Lhat_RD(delta, Cj, Cbar, Xt, Xc, mon_ind,
                               sigma_t, sigma_c, ht, hc)

  sd <- sd_Lhat_RD(delta, Cj, Cbar, Xt, Xc, mon_ind,
                   sigma_t, sigma_c, ht, hc)

  res <- list(ci.l = lhat - sup_bias - stats::qnorm(1 - tau) * sd,
              sd = sd)
  return(res)
}
