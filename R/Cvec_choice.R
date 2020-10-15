#' Calculate \eqn{E U_j} in the paper
#'
#' @param Cpr \eqn{C'} in the paper
#' @param C the Lipschitz coefficient for the function space we consider
#' @param tau_res a list produced by \code{tau_calc}; can be left unspecified
#' @param hmat a matrix of bandwidths to be used in the adaptive procedure;
#' can be left unspecified
#' @inheritParams tau_calc
#'
#' @return J-dimensional vector for \eqn{E U_j}'s
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

  if(missing(tau_res)){
    tau_res <- tau_calc(Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha)
  }

  tau_sol <- tau_res$tau_sol
  del_sol <- tau_res$del_sol

  J <- length(Cvec)

  Yt <- -Cpr * Norm(Vminus(Xt, mon_ind))
  Yc <- Cpr * Norm(Vplus(Xc, mon_ind))

  EU_vec <- numeric(J)

  for(j in 1:J){

    Cj <- Cvec[j]

    if(missing(hmat)){

      hres_j <- bw_adpt(del_sol, Cj, C, Xt, Xc, mon_ind, sigma_t, sigma_c)
      ht_j <- hres_j$ht
      hc_j <- hres_j$hc

    }else{

      ht_j <- hmat[j, 1]
      hc_j <- hmat[j, 2]
    }

    EU_vec[j] <- -c_hat_lower_RD(del_sol, Cj, C, Xt, Xc, mon_ind,
                             sigma_t, sigma_c, Yt, Yc, tau_sol, ht_j, hc_j)
  }

  return(EU_vec)
}


#' Calculate \eqn{Cov(U_j, U_k)} in the paper
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
      ht_j <- hres_j$ht
      hc_j <- hres_j$hc

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

#' Row-wise minimum
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

#' Worst-case Excess Length of the Adaptive Procedure
#'
#' @inheritParams EU_vec
#' @param sim number of simulated observations to calculate the
#' expectation of the minimum of multivariate normal random variables.
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
                   sim = 10^5, tau_res, hmat){

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
      hmat[j, 1] <- hres_j$ht
      hmat[j, 2] <- hres_j$hc
    }
  }

  EUs <- EU_vec(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                tau_res, hmat)
  CovU <- CovU_mat(Cpr, Cvec, C, Xt, Xc, mon_ind, sigma_t, sigma_c, alpha,
                   tau_res, hmat)

  simvec <- mvtnorm::rmvnorm(sim, EUs, CovU)
  minvec <- pmin_mat(simvec)
  res <- mean(minvec)

  return(res)
}

#' Oracle Excess Length
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

  C_pair <- c(C, Cpr)
  modres <- modsol_RD(stats::qnorm(1 - alpha), C_pair, Xt, Xc, mon_ind,
                      sigma_t, sigma_c)

  res <- modres$bt + modres$bc
  return(res)
}
