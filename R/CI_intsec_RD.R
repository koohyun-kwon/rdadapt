#' Covariance Between Two Estimators
#'
#' Calculates the covariance between two estimators,
#' \eqn{Lhat_j(\delta)} and \eqn{Lhat_k(\delta)}.
#'
#' This corresponds to the expressions (17) and (18) of our paper.
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
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' cov_calc(1, 0.2, 0.4, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' cov_calc(1, 0.2, 0.4, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
cov_calc <- function(delta, Cj, Ck, Cbar, Xt, Xc, mon_ind,
                     sigma_t, sigma_c, ht_j, hc_j, ht_k, hc_k){

  if(missing(ht_j) | missing(hc_j) | missing(ht_k) | missing(hc_k)){

    bres_j <- bw_mod(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_j <- bres_j$bt
    hc_j <- bres_j$bc

    bres_k <- bw_mod(delta, Ck, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
    ht_k <- bres_k$bt
    hc_k <- bres_k$bc
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
#' This constructs a covariance matrix corresponding to
#' expressions (17) and (18) of our paper.
#'
#' @param delta a nonegative scalar value;
#' it can be left unspecified if \code{bmat} is specified.
#' @param Cvec a sequence of smoothness parameters
#' @param Cbar the Lipschitz coefficient for the largest function space we consider
#' @inheritParams cov_calc
#' @param bmat a \eqn{J} by {2} matrix of modulus values;
#' it can be left unspecified if \code{delta} is specified.
#'
#' @return a \eqn{J} by \eqn{J} covariance matrix.
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
#' cov_mat_calc(1, (1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c)
#' cov_mat_calc(1, (1:5)/5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c)
cov_mat_calc <- function(delta, Cvec, Cbar, Xt, Xc, mon_ind,
                         sigma_t, sigma_c, bmat){

  J <- length(Cvec)

  if(missing(bmat)){

    bmat <- matrix(0, nrow = J, ncol = 2)

    for(j in 1:J){

      Cj <- Cvec[j]

      bres_j <- bw_mod(delta, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      bmat[j, 1] <- bres_j$bt
      bmat[j, 2] <- bres_j$bc
    }
  }

  res <- diag(1, J, J)

  if(J > 1){

    for(j in 1:(J - 1)){

      for(k in (j + 1):J){

        Cj <- Cvec[j]
        Ck <- Cvec[k]
        ht_j <- bmat[j, 1]
        hc_j <- bmat[j, 2]
        ht_k <- bmat[k, 1]
        hc_k <- bmat[k, 2]

        res[j, k] <- cov_calc(delta, Cj, Ck, Cbar, Xt, Xc, mon_ind,
                              sigma_t, sigma_c, ht_j, hc_j, ht_k, hc_k)
        res[k, j] <- res[j, k]
      }
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
#' the default is \code{10^5}.
#'
#' @return \eqn{(1 - \alpha)}th quantile of maximum of a multivariate normal random variable.
#' @export
#'
#' @examples covmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' max_Q(covmat, 0.05)
max_Q <- function(covmat, alpha, num_sim = 10^5){

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
#' This function was created for a test purpose; it is not used since it is slow for
#' large \eqn{J}.
#'
#' @param covmat covariance matrix of the multivariate normal random variable.
#' @param alpha desired upper quantile value.
#'
#' @return \eqn{(1 - \alpha)}th quantile of maximum of a multivariate normal random variable.
#' @export
#'
#' @examples covmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' max_Q2(covmat, 0.05)
max_Q2 <- function(covmat, alpha){

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
#' This solves (19) of our paper, using the asymptotic argument provided in the paper.
#'
#' @inheritParams cov_mat_calc
#' @inheritParams max_Q
#' @param delta_init the value of \eqn{\delta} to be used in simulating the quantile;
#' theoretically, its value does not matter asymptotically. Its default value is 1.96.
#' @param bmat_init the matrix of modulus values corresponding to \code{delta_init}
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
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' tau_calc((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
#' tau_calc((1:5)/5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
tau_calc <- function(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                    num_sim = 10^5, delta_init = 1.96, bmat_init){

  J <- length(Cvec)

  if(missing(bmat_init)){

    bmat_init <- matrix(0, nrow = J, ncol = 2)

    for(j in 1:J){

      Cj <- Cvec[j]

      bres_j <- bw_mod(delta_init, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      bmat_init[j, 1] <- bres_j$bt
      bmat_init[j, 2] <- bres_j$bc
    }
  }

  covmat <- cov_mat_calc(delta_init, Cvec, Cbar, Xt, Xc, mon_ind,
                         sigma_t, sigma_c, bmat_init)
  del_sol <- max_Q(covmat, alpha, num_sim)
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
#' @param bmat the matrix of modulus values corresponding to
#' the optimal \eqn{\delta} and \code{Cvec}; it can be left unspecified.
#' @param tau_res outcome value from the function \code{tau_res}; it can be
#' left unspecified.
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
#' tind <- X[, 1] < 0 & X[, 2] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' CI_adpt_L((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05)
#' CI_adpt_L((1:5)/5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05)
CI_adpt_L <- function(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                      num_sim = 10^5, delta_init = 1.96, bmat_init, bmat, tau_res){

  J <- length(Cvec)

  if(missing(bmat_init)){

    if(missing(bmat)){

      bmat_init <- matrix(0, nrow = J, ncol = 2)

      for(j in 1:J){

        Cj <- Cvec[j]

        bres_j <- bw_mod(delta_init, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
        bmat_init[j, 1] <- bres_j$bt
        bmat_init[j, 2] <- bres_j$bc
      }
    }
  }

  if(missing(tau_res)){

    tau_res <- tau_calc(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha, num_sim,
                        delta_init, bmat_init)
  }

  tau_sol <- tau_res$tau_sol
  del_sol <- tau_res$del_sol

  resvec <- numeric(J)

  for(j in 1:J){

    Cj <- Cvec[j]

    if(missing(bmat)){

      bres_j <- bw_mod(del_sol, Cj, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c)
      ht_j <- bres_j$bt
      hc_j <- bres_j$bc

    }else{

      ht_j <- bmat[j, 1]
      hc_j <- bmat[j, 2]

    }

    resvec[j] <- c_hat_lower_RD(del_sol, Cj, Cbar, Xt, Xc, mon_ind,
                                sigma_t, sigma_c, Yt, Yc, tau_sol,
                                ht_j, hc_j)$ci.l
  }

  res_list <- list(CI = max(resvec), CI_vec = resvec, tau_sol = tau_sol, del_sol = del_sol,
                   Cj = Cvec[which.max(resvec)])
  return(res_list)
}


#' Adaptive Confidence Interval
#'
#' Calculates the one-sided or two-sided adaptive confidence interval
#'
#' This returns either one-sided lower CI or two-sided Bonferroni CI using
#' the function CI_adpt_L; it seems this function will be substituted by
#' \code{CI_adpt_opt} defined below.
#'
#' @inheritParams CI_adpt_L
#' @param lower calculate a lower one-sided confidence interval if \code{TRUE};
#' calculate a two-sided CI otherwise.
#' @param bmat_init_L the matrix of modulus values corresponding to
#' \code{delta_init} and \code{Cvec}; it can be left unspecified.
#' @param bmat_L the matrix of modulus values corresponding to
#' the optimal \eqn{\delta} and \code{Cvec}; it can be left unspecified.
#'
#' @return confidence interval endpoints.
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
#' CI_adpt((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = FALSE)
#' CI_adpt((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = TRUE)
#' CI_adpt((1:5)/5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = TRUE)
CI_adpt <- function(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                    lower = TRUE, num_sim = 10^5, delta_init = 1.96, bmat_init_L,
                    bmat_L, tau_res){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  if(lower == T){

    res_L <- CI_adpt_L(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, num_sim,
                     delta_init = 1.96, bmat_init_L, bmat_L, tau_res)
    res <- c(res_L$CI, Inf)
  }else{

    res_L <- CI_adpt_L(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha/2, num_sim,
                       delta_init = 1.96, bmat_init_L, bmat_L, tau_res)
    CI_L <- res_L$CI

    CI_U <- -c_hat_lower_RD(stats::qnorm(1 - alpha/2), Cbar, Cbar, Xc, Xt, mon_ind,
                            sigma_c, sigma_t, Yc, Yt, alpha/2)$ci.l

    res <- c(CI_L, CI_U)
  }

  return(res)
}


#' Adaptive Confidence Interval with Optimal Lipschitz Coefficients
#'
#' Constructs an adaptive lower CI after choosing the optimal sequence of
#' Lipschitz coefficients.
#'
#' This function also supports variance estimation.
#'
#' @inheritParams Opt_C_seq
#' @inheritParams CI_adpt
#' @param sigma_t.init supplied first-stage variance for treated observations.
#' @param sigma_c.init supplied first-stage variance for control observations.
#' @param se.method standard deviation estimation methods.
#' @param J a positive integer; if specified,
#' the sequence of Lipschitz coefficients is set to be J equidistant grids
#' in \code{(C_l, C_u)}, without the optimal choice procedure.
#' @param se.init the standard deviation estimation method for choosing an optimal estimator.
#' @param t.dir treatment direction; \code{t.dir = "left"} if \eqn{x < 0} is treated.
#' Otherwise, \code{t.dir = "right"}. This should specified only for one-dimensional cases.
#' @param N number of nearest neighbors to match when doing the variance estimation.
#'
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
#' CI_adpt_opt(0.1, 1, 2, Xt, Xc, mon_ind, Yt, Yc, 0.05, "supplied", sigma_t, sigma_c, J = 5,
#' se.init = "supplied")
#' d <- 1
#' X <- rnorm(n)
#' tind <- X < 0
#' Xt <- X[tind == 1]
#' Xc <- X[tind == 0]
#' mon_ind <- 1
#' sigma <- rep(1, n)
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt <- 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' CI_adpt_opt(0.1, 1, 2, Xt, Xc, mon_ind, Yt, Yc, 0.05, "nn", se.init = "Silverman",
#' t.dir = "left")
#' CI_adpt_opt(1, 1, 2, Xt, Xc, mon_ind, Yt, Yc, 0.05, "nn", se.init = "Silverman",
#' t.dir = "left")
CI_adpt_opt <- function(C_l, C_u, C, Xt, Xc, mon_ind, Yt, Yc, alpha,
                        se.method = c("nn", "supplied", "nn.test"), sigma_t, sigma_c,
                        sigma_t.init, sigma_c.init, J,
                        se.init = c("Silverman", "supplied", "supp.sep", "S.test"),
                        t.dir = c("left", "right"),
                        n_grid = 10,
                        gain_tol = 0.05, ratio = TRUE, p = Inf, n_sim = 10^5,
                        delta_init = 1.96, N = 3){

  # Convert running variable data into matrices
  if(!is.matrix(Xt)) Xt <- matrix(Xt, ncol = 1)
  if(!is.matrix(Xc)) Xc <- matrix(Xc, ncol = 1)

  se.method = match.arg(se.method)
  se.init = match.arg(se.init)
  t.dir = match.arg(t.dir)

  # The sorting process is necessary due to the way variance estimation functions work
  # Only supports d = 1 case
  if(ncol(Xt) == 1){
    sorted.t <- sort(Xt, index.return = T)
    sorted.c <- sort(Xc, index.return = T)

    Xt <- matrix(sorted.t$x, ncol = 1)
    Xc <- matrix(sorted.c$x, ncol = 1)
    Yt <- Yt[sorted.t$ix]
    Yc <- Yc[sorted.c$ix]

    if(!missing(sigma_t) & !missing(sigma_c)){

      sigma_t <- sigma_t[sorted.t$ix]
      sigma_c <- sigma_c[sorted.c$ix]
    }

    if(!missing(sigma_t.init) & !missing(sigma_c.init)){

      sigma_t.init <- sigma_t.init[sorted.t$ix]
      sigma_c.init <- sigma_c.init[sorted.c$ix]
    }
  }

  if(se.init == "supplied"){

    sigma_t.init <- sigma_t
    sigma_c.init <- sigma_c

  }else if(se.init == "Silverman"){

    if(ncol(Xt) > 1){
      stop("Multi-dimension not supported for now")
    }

    sigma.init <- sigmaSvm(Xt, Xc, Yt, Yc, t.dir)
    sigma_t.init <- sigma.init$sigma.t
    sigma_c.init <- sigma.init$sigma.c
  }else if(se.init == "S.test"){

    sigma.init <- sigmaSvm.test(Xt, Xc, Yt, Yc)
    sigma_t.init <- sigma.init$sigma.t
    sigma_c.init <- sigma.init$sigma.c
  }

  if(missing(J) & C_l != C_u){

    opt_Cvec_res <- Opt_C_seq(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t.init, sigma_c.init, alpha,
                              n_grid, gain_tol, ratio, p, n_sim)

    Cvec <- opt_Cvec_res$Cvec
    bmat <- opt_Cvec_res$bmat
    tau_res <- opt_Cvec_res$tau_res

    res_L <- CI_adpt_L(Cvec, C, Xt, Xc, mon_ind, sigma_t.init, sigma_c.init, Yt, Yc,
                       alpha, n_sim, delta_init, bmat = bmat, tau_res = tau_res)

  }else{

    if(C_l == C_u) J = 1

    C_grid <- seq(from = C_l, to = C_u, length.out = J + 2)
    Cvec <- C_grid[2:(J + 1)]

    res_L <- CI_adpt_L(Cvec, C, Xt, Xc, mon_ind, sigma_t.init, sigma_c.init, Yt, Yc,
                       alpha, n_sim, delta_init)
  }


  # sd estimation used to calculate the final CI
  # Note that  if se.method == "supplied", sigma_t, sigma_c values are given
  if(se.method == "nn"){

    if(ncol(Xt) > 1){
      stop("Multi-dimension not supported for now")
    }

    sigma.res <- sigmaNN(Xt, Xc, Yt, Yc, t.dir, N)
    sigma_t.new <- sigma.res$sigma.t
    sigma_c.new <- sigma.res$sigma.c

  }else if(se.method == "nn.test"){

    sigma.res <- sigmaNN2(Xt, Xc, Yt, Yc, N)
    sigma_t.new <- sigma.res$sigma.t
    sigma_c.new <- sigma.res$sigma.c

  }else if(se.method == "supplied"){

    sigma_t.new <- sigma_t
    sigma_c.new <- sigma_c
  }

  C.chosen <- res_L$Cj

  tau_sol <- res_L$tau_sol
  del_sol <- res_L$del_sol

  # Modifies the CI so that it uses the 2nd stage sd estimator
  sd.old <- c_hat_lower_RD(del_sol, C.chosen, C_u, Xt, Xc, mon_ind,
                          sigma_t.init, sigma_c.init, Yt, Yc, tau_sol)$sd

  lhat.w_t <- Lhat_fun_RD(del_sol, C.chosen, C_u, Xt, Xc, mon_ind,
                        sigma_t.init, sigma_c.init, Yt, Yc, ret.w = TRUE)$w_t
  lhat.w_c <- Lhat_fun_RD(del_sol, C.chosen, C_u, Xt, Xc, mon_ind,
                          sigma_t.init, sigma_c.init, Yt, Yc, ret.w = TRUE)$w_c
  sd.new <- sd_w(lhat.w_t, lhat.w_c, sigma_t.new, sigma_c.new)

  res.l <- res_L$CI - sd.old * stats::qnorm(1 - tau_sol) + sd.new * stats::qnorm(1 - tau_sol)

  res <- c(res.l, Inf)
  return(res)
}
