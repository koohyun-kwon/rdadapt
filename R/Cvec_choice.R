#' Distance from Oracle Helper Function: Mean
#'
#' Calculates \eqn{E U_j} in Proposition 4.4 of the paper (Ver20201104)
#'
#' @param Cpr the Lipschitz coefficient where the distance is evaluated at.
#' @param C the Lipschitz coefficient for the function space we consider.
#' @param tau_res a list produced by the function \code{tau_calc};
#' can be left unspecified.
#' @param hmat a matrix of bandwidths to be used in the adaptive procedure;
#' can be left unspecified.
#' @inheritParams tau_calc
#'
#' @return J-dimensional vector containing values of \eqn{E U_j}'s.
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
#' EU_vec(1/4, (1:5)/5, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
EU_vec <- function(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                   tau_res, hmat){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  if(missing(tau_res)){
    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  tau_sol <- tau_res$tau_sol
  del_sol <- tau_res$del_sol

  J <- length(Cvec)

  # "Pseudo-" Yt and Yc values
  pYt <- -Cpr * Norm(Vminus(Xt, mon_ind))
  pYc <- Cpr * Norm(Vplus(Xc, mon_ind))

  EU_vec <- numeric(J)

  for(j in 1:J){

    Cj <- Cvec[j]

    if(missing(hmat)){

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      ht_j <- hres_j$bt
      hc_j <- hres_j$bc

    }else{

      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]
    }

    EU_vec[j] <- -c_hat_lower_RD(del_sol, Cj, C, Xt, Xc, mon_ind,
                             sigma_t, sigma_c, pYt, pYc, tau_sol, ht_j, hc_j)$ci.l
  }

  return(EU_vec)
}


#' Distance from Oracle Helper Function: Covariance
#'
#' Calculates \eqn{Cov(U_j, U_k)} in Proposition 4.4 of the paper (Ver20201104).
#'
#' @inheritParams EU_vec
#'
#' @return Covariance matrix with dimension J by J
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
#' CovU_mat(1/4, (1:5)/5, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
CovU_mat <- function(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                     tau_res, hmat){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  if(missing(tau_res)){
    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  tau_sol <- tau_res$tau_sol
  del_sol <- tau_res$del_sol

  J <- length(Cvec)
  nt <- nrow(Xt)
  nc <- nrow(Xc)

  w_mat_t <- matrix(0, nrow = J, ncol = nt)
  w_mat_c <- matrix(0, nrow = J, ncol = nc)

  for(j in 1:J){

    Cj <- Cvec[j]

    if(missing(hmat)){

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      ht_j <- hres_j$bt
      hc_j <- hres_j$bc

    }else{

      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]
    }

    num_tj <- K_fun(b = ht_j, C_pair = c(C, Cj), Xt, mon_ind) / sigma_t
    num_cj <- K_fun(b = hc_j, C_pair = c(Cj, C), Xc, mon_ind) / sigma_c

    denom_tj <- sum(K_fun(b = ht_j, C_pair = c(C, Cj), Xt, mon_ind) / sigma_t^2)
    denom_cj <- sum(K_fun(b = hc_j, C_pair = c(Cj, C), Xc, mon_ind) / sigma_c^2)

    w_mat_t[j, ] <- num_tj / denom_tj
    w_mat_c[j, ] <- num_cj / denom_cj
  }

  res <- w_mat_t %*% t(w_mat_t) + w_mat_c %*% t(w_mat_c)
  return(res)
}

#' Worst-case Excess Length of the Adaptive Procedure
#'
#' Calculates E[U_{min}] in Proposition 4.4 of the paper (Ver20201104).
#'
#' @inheritParams EU_vec
#' @param n_sim number of simulated observations to calculate the
#' expectation of the minimum of multivariate normal random variables;
#' the default is \code{n_sim = 10^5}.
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
#' l_adpt(1/4, (1:5)/5, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
l_adpt <- function(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                   n_sim = 10^5, tau_res, hmat){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  J <- length(Cvec)

  if(missing(tau_res)){
    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  if(missing(hmat)){

    del_sol <- tau_res$del_sol

    hmat <- matrix(0, nrow = J, ncol = 2)
    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$bt
      hmat[j, 2] <- hres_j$bc
    }
  }

  EUs <- EU_vec(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                tau_res, hmat)
  CovU <- CovU_mat(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                   tau_res, hmat)

  simvec <- mvtnorm::rmvnorm(n_sim, EUs, CovU)
  minvec <- pmin_mat(simvec)
  res <- mean(minvec)

  return(res)
}

#' Oracle Excess Length
#'
#' Calculates the expression for \eqn{\ell(C';C)}
#' in Proposition 4.4 of the paper (Ver20201104).
#'
#' @inheritParams l_adpt
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
#' l_orc(1/4, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
l_orc <- function(Cpr, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  C_pair <- c(C, Cpr)
  modres <- modsol_RD(stats::qnorm(1 - alpha), C_pair, Xt, Xc, mon_ind,
                      sigma_t, sigma_c)

  res <- modres$bt + modres$bc
  return(res)
}

#' Distance Between Two Length Functions
#'
#' Calculates \eqn{\Delta(l_1, l_2)} in the paper (Ver20201104).
#'
#' @param l_adpt_vec vector of worst case lengths for the adaptive procedure.
#' @param l_orc_vec vector of worst case lengths for the oracle procedure.
#' @param ratio the ratio measure is used if \code{TRUE};
#' otherwise, the difference measure is used.
#' @param p the order of \eqn{l_1}-norm; the default is \code{Inf}.
#'
#' @export
#'
#' @examples l_adpt_vec <- 1:10
#' l_orc_vec <- (1:10) - 1/2
#' dist_l(l_adpt_vec, l_orc_vec)
dist_l <- function(l_adpt_vec, l_orc_vec, ratio = TRUE, p = Inf){

  if(ratio == TRUE){

    l_diff <- l_adpt_vec / l_orc_vec

  }else{

    l_diff <- l_adpt_vec - l_orc_vec
  }

  if(p == Inf){

    res <- max(l_diff)

  }else{

    res <- (sum(l_diff^p))^(1 / p)
  }

  return(res)
}

#' Calculates the Optimal Sequence of Lipschitz Coefficients
#'
#' This corresponds to Procedure 3 in the paper (Ver20201104).
#'
#' @param C_l lower end of the adaptation range.
#' @param C_u upper end of the adaptation range.
#' @param n_grid number of grid points to evaluate the lengths.
#' @param gain_tol stopping criterion when finding the optimal J.
#' @inheritParams l_adpt
#' @inheritParams dist_l
#'
#' @return a list with components \code{Cvec}, the optimal sequence of Lipschitz coefficients,
#' \code{hmat}, the matrix of corresponding bandwidths,
#' \code{tau_res}, the corresponding calibrated values of \eqn{\tau} and \eqn{\delta}, and
#' \code{dist_opt}, the optimal distance to the orale.
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
#' Opt_C_seq(0.1, 1, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05, 10, 0.05)
Opt_C_seq <- function(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                      n_grid = 10, gain_tol = 0.05, ratio = TRUE, p = Inf, n_sim = 10^5){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  # The length measure is evaluated over this grid
  Cpr_grid <- seq(from = C_l, to = C_u, length.out = n_grid)

  # Calculates worst-case length of the oracle procedure over Cpr_grid
  l_orc_vec <- numeric(n_grid)
  for(i in 1:n_grid){

    l_orc_vec[i] <- l_orc(Cpr_grid[i], C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  J <- 0
  dist <- Inf # initialization of the distance result
  gain <- 1 # initialization of the gain from increasing J

  # The loop will stop when relative distance reduction from adapting to one more space
  # is less than gain_tol
  while(gain > gain_tol){

    J <- J + 1
    C_grid <- seq(from = C_l, to = C_u, length.out = J + 2)
    Cvec <- C_grid[2:(J + 1)] # adapt to J middle points over [C_l, C_u]

    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
    del_sol <- tau_res$del_sol

    hmat <- matrix(0, nrow = J, ncol = 2)
    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$bt
      hmat[j, 2] <- hres_j$bc
    }

    # Calculates worst-case length of the adaptive procedure over Cpr_grid
    l_adpt_vec <- numeric(n_grid)
    for(i in 1:n_grid){

      l_adpt_vec[i] <- l_adpt(Cpr_grid[i], Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                              alpha, n_sim, tau_res, hmat)
    }

    dist_new <- dist_l(l_adpt_vec, l_orc_vec, ratio, p)
    gain <- 1 - dist_new / dist # percentage gain
    dist <- dist_new
  }

  res <- list(Cvec_opt = Cvec, hmat = hmat, tau_res = tau_res,
              dist_opt = dist)

  return(res)

}


#' Distance Measure Function (of J)
#'
#' @inheritParams Opt_C_seq
#' @param J the value of J
#'
#' @return the distance measure
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
#' dist.J(0.1, 1, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05, 2)
dist.J <- function(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha, J,
                   n_grid = 10, ratio = TRUE, p = Inf, n_sim = 10^5){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  # The length measure is evaluated over this grid
  Cpr_grid <- seq(from = C_l, to = C_u, length.out = n_grid)

  # Calculates worst-case length of the oracle procedure over Cpr_grid
  l_orc_vec <- numeric(n_grid)
  for(i in 1:n_grid){

    l_orc_vec[i] <- l_orc(Cpr_grid[i], C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  C_grid <- seq(from = C_l, to = C_u, length.out = J + 2)
  Cvec <- C_grid[2:(J + 1)] # adapt to J middle points over [C_l, C_u]

  tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  del_sol <- tau_res$del_sol

  hmat <- matrix(0, nrow = J, ncol = 2)
  for(j in 1:J){

    Cj <- Cvec[j]

    hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
    hmat[j, 1] <- hres_j$bt
    hmat[j, 2] <- hres_j$bc
  }

  # Calculates worst-case length of the adaptive procedure over Cpr_grid
  l_adpt_vec <- numeric(n_grid)
  for(i in 1:n_grid){

    l_adpt_vec[i] <- l_adpt(Cpr_grid[i], Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                            alpha, n_sim, tau_res, hmat)
  }

  dist <- dist_l(l_adpt_vec, l_orc_vec, ratio, p)

  return(dist)
}



#' Calculates the Optimal Sequence of Lipschitz Coefficients 2
#'
#' This corresponds to Procedure 3 in the paper (Ver20201104).
#'
#' @param C_l lower end of the adaptation range.
#' @param C_u upper end of the adaptation range.
#' @param gain_tol stopping criterion when finding the optimal J.
#' @inheritParams l_adpt
#' @inheritParams dist_l
#'
#' @return a list with components \code{Cvec}, the optimal sequence of Lipschitz coefficients,
#' \code{hmat}, the matrix of corresponding bandwidths,
#' \code{tau_res}, the corresponding calibrated values of \eqn{\tau} and \eqn{\delta}, and
#' \code{dist_opt}, the optimal distance to the orale.
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
#' Opt_C_seq2(0.1, 1, 2, Xt, Xc, mon_ind, sigma_t, sigma_c, 0.05)
Opt_C_seq2 <- function(C_l, C_u, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                      gain_tol = 0.05, ratio = TRUE, p = Inf, n_sim = 10^5){

  if(!is.matrix(Xt) | !is.matrix(Xc)){
    stop("Xt and Xc should be matrices")
  }

  J <- 0
  dist <- Inf # initialization of the distance result
  gain <- 1 # initialization of the gain from increasing J

  # The loop will stop when relative distance reduction from adapting to one more space
  # is less than gain_tol
  while(gain > gain_tol){

    J <- J + 1
    C_grid <- seq(from = C_l, to = C_u, length.out = J + 2)
    Cvec <- C_grid[2:(J + 1)] # adapt to J middle points over [C_l, C_u]

    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
    del_sol <- tau_res$del_sol

    hmat <- matrix(0, nrow = J, ncol = 2)
    for(j in 1:J){

      Cj <- Cvec[j]

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      hmat[j, 1] <- hres_j$bt
      hmat[j, 2] <- hres_j$bc
    }

    # Calculates the worst-case distance
    dist_fun <- function(Cpr){

      l_adpt(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
             alpha, n_sim, tau_res, hmat) /
        l_orc(Cpr, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
    }

    dist.res <- stats::optimize(dist_fun, c(C_l, C_u), maximum = TRUE)

    dist_new <- dist.res$objective
    gain <- 1 - dist_new / dist # percentage gain
    dist <- dist_new
  }

  res <- list(Cvec_opt = Cvec, hmat = hmat, tau_res = tau_res,
              dist_opt = dist)

  return(res)

}
