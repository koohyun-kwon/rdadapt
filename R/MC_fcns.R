#' Monte Carlo Simulation for Adaptive Lower CI
#'
#' @inheritParams Opt_C_seq
#' @inheritParams CI_adpt_opt
#' @param Yt_mat a \eqn{n_t} by \eqn{m} outcome matrix, where \eqn{m} is the number of simulated
#' observations.
#' @param Yc_mat  a \eqn{n_c} by \eqn{m} outcome matrix, where \eqn{m} is the number of simulated
#' observations.
#' @param val_tr the true value of the parameter
#'
#' @return a vector of the coverage probability and the average length.
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
#' val_tr <- 1
#' Yt1 <- val_tr + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yt2 <- val_tr + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc1 <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' Yc2 <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' Yt_mat <- cbind(Yt1, Yt2)
#' Yc_mat <- cbind(Yc1, Yc2)
#' MC_sim_lower(0.1, 1, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt_mat, Yc_mat, val_tr,
#' 0.05)
#' MC_sim_lower(0.1, 1, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt_mat, Yc_mat, val_tr,
#' 0.05, 3)
MC_sim_lower <- function(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt_mat, Yc_mat, val_tr,
                   alpha, J, n_grid = 10, gain_tol = 0.05, ratio = TRUE, p = Inf, n_sim = 10^5,
                   delta_init = 1.96){

  if(missing(J)){

    opt_Cvec_res <- Opt_C_seq(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                              n_grid, gain_tol, ratio, p, n_sim)

    Cvec <- opt_Cvec_res$Cvec
    hmat <- opt_Cvec_res$hmat
    tau_res <- opt_Cvec_res$tau_res

  }else{

    C_grid <- seq(from = C_l, to = C_u, length.out = J + 2)
    Cvec <- C_grid[2:(J + 1)]

    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
    del_sol <- tau_res$del_sol

    hmat <- matrix(0, nrow = J, ncol = 2)
    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$ht
      hmat[j, 2] <- hres_j$hc
    }
  }

  n_MC <- ncol(Yt_mat)

  cov_res <- numeric(n_MC)
  len_res <- numeric(n_MC)

  for(i in 1:n_MC){

    Yt <- Yt_mat[, i]
    Yc <- Yc_mat[, i]

    res_L <- CI_adpt_L(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, n_sim,
                       hmat = hmat, tau_res = tau_res)

    cov_res[i] <- res_L$CI[1] <= val_tr
    len_res[i] <- val_tr - res_L$CI[1]
  }

  res <- c(mean(cov_res), mean(len_res))
  names(res) <- c("cov_prob", "avg_len")

  return(res)
}

#' Monte Carlo Simulations for Minimax Procedure
#'
#' @inheritParams MC_sim_lower
#' @inheritParams CI_minimax_RD
#'
#' @return a vector of the coverage probability and the minimax length.
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
#' val_tr <- 1
#' Yt1 <- val_tr + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yt2 <- val_tr + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc1 <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' Yc2 <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' Yt_mat <- cbind(Yt1, Yt2)
#' Yc_mat <- cbind(Yc1, Yc2)
#' MC_sim_mm(Yt_mat, Yc_mat, Xt, Xc, 2, mon_ind, sigma_t, sigma_c, 0.05, val_tr)
MC_sim_mm <- function(Yt_mat, Yc_mat, Xt, Xc, C, mon_ind, sigma_t, sigma_c, alpha, val_tr,
                      maxb.const = 10){

  minbt <- minb_fun(rep(C, 2), Xt, mon_ind)
  minbc <- minb_fun(rep(C, 2), Xc, mon_ind, swap = T)
  minb <- minbt + minbc

  modres_2 <- modsol_RD(stats::qnorm(1 - alpha/2), rep(C,2), Xt, Xc,
                        mon_ind, sigma_t, sigma_c)

  maxb <- maxb.const * (modres_2$bt + modres_2$bc)

  CI_length_sol <- stats::optimize(CI_length_RD, interval = c(minb, maxb),
                                   C = C, Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                                   sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)

  min_half_length <- CI_length_sol$objective
  opt_b <- CI_length_sol$minimum

  invmod_opt <- invmod_RD(opt_b, rep(C, 2), Xt, Xc, mon_ind, sigma_t, sigma_c)
  ht <- invmod_opt$bt
  hc <- invmod_opt$bc
  del_t <- invmod_opt$delta_t
  del_c <- invmod_opt$delta_c
  delta <- sqrt(del_t^2 + del_c^2)

  a_val <- a_fun(delta, C, C, Xt, Xc, mon_ind, sigma_t, sigma_c, ht, hc)

  n_MC <- ncol(Yt_mat)

  cov_res <- numeric(n_MC)

  for(i in 1:n_MC){

    Yt <- Yt_mat[, i]
    Yc <- Yc_mat[, i]
    opt_Lhat <- Lhat_fun_RD(delta, C, C, Xt, Xc, mon_ind,
                            sigma_t, sigma_c, Yt, Yc, ht, hc) - a_val

    CIres <- c(opt_Lhat - min_half_length, opt_Lhat + min_half_length)

    cov_res[i] <- (val_tr >= CIres[1]) & (val_tr <= CIres[2])
  }

  res <- c(mean(cov_res), 2 * min_half_length)
  names(res) <- c("cov_prob", "avg_len")

  return(res)
}
