#' Specification Setting for Simulation
#'
#' @param spec_num a number symbolizing the specification to be used.
#' @param trueC_len the number of true smoothness spaces.
#' @param trueC_sc a scalar value of the smoothness parameter,
#' used for power function calculations.
#'
#'
#' @return a list of elements.
#' @export
#'
#' @examples spec_set(1, 5)
#' spec_set(2, 5)
#' spec_set(1, 1, 0.5)
spec_set <- function(spec_num, trueC_len, trueC_sc = NULL){

  if(spec_num == 1){  # Linear, 1-dim

    trueC_min <- 1/5
    trueC_max <- 1
    trueC <- seq(from = trueC_min, to = trueC_max, length.out = trueC_len)

    d <- 1
    mon_ind <- c(1)

    spec <- "Lin"
    sd_spec <- "hom"
    X_dist = "equi"

    alpha <- 0.05

  }else if(spec_num == 2){ # Linear, 2-dim

    trueC_min <- 1/5
    trueC_max <- 1
    trueC <- seq(from = trueC_min, to = trueC_max, length.out = trueC_len)

    d <- 2
    mon_ind <- c(1, 2)

    spec <- "Lin"
    sd_spec <- "hom"
    X_dist = "equi"

    alpha <- 0.05

  }else if(spec_num == 3){ # Linear, 1-dim, larger C_max

    trueC_min <- 1/5
    trueC_max <- 5
    trueC <- seq(from = trueC_min, to = trueC_max, length.out = trueC_len)

    d <- 1
    mon_ind <- c(1)

    spec <- "Lin"
    sd_spec <- "hom"
    X_dist = "equi"

    alpha <- 0.05
  }

  return(list(d = d, mon_ind = mon_ind, trueC = trueC, fcn_spec = spec, sd_spec = sd_spec,
              X_dist = X_dist, alpha = alpha, trueC_sc = trueC_sc))
}


#' Sample Generation for Simulations
#'
#' @param n sample size.
#' @param d number of running variables.
#' @param C a scalar or a vector of true smoothness parameters.
#' @param spec specifications to use.
#' @param sd_spec homoskedastic error if \code{sd_spec = "hom"};
#' heteroskedastic if \code{sd_spec = "het"}.
#' @param X_dist distribution of running variables
#' @param true_val the true parameter value; the default is 1.
#'
#' @return a list with \code{Yt}, \code{Yc}, \code{Xt}, \code{Xc}, \code{sigma_t},
#' and \code{sigma_c}.
#' @export
#'
#' @examples gen_obs(100, 1, 1/2, "Lin", "het", "equi")
#' gen_obs(100, 1, c(1/2, 1), "Lin", "hom", "equi")
#' gen_obs(200, 2, c(1/2, 1), "Lin", "hom", "equi")
gen_obs <- function(n, d, C, spec = c("Lin"), sd_spec = c("hom", "het"),
                    X_dist = c("equi"), true_val = 1){

  X_dist <- match.arg(X_dist)
  sd_spec <- match.arg(sd_spec)
  spec <- match.arg(spec)
  trueC_len <- length(C)

  if(sd_spec == "hom"){

    sig = rep(1, n)

  }else if(sd_spec == "het"){

    sig = rep(1, n) #to be revised
  }

  if(d == 1){

    if(X_dist == "equi"){

      X <- matrix(seq(from = -1, to = 1, length.out = n), nrow = n, ncol = d)
      tind <- X < 0
    }

    Xt <- X[tind == 1, ,drop = F]
    Xc <- X[tind == 0, ,drop = F]
    sigma_t <- sig[tind == 1]
    sigma_c <- sig[tind == 0]

    if(spec == "Lin"){

      ut <- rep(stats::rnorm(length(sigma_t), mean = 0, sd = sigma_t), trueC_len)
      ut <- matrix(ut, nrow = length(sigma_t), ncol = trueC_len)
      uc <- rep(stats::rnorm(length(sigma_c), mean = 0, sd = sigma_c), trueC_len)
      uc <- matrix(uc, nrow = length(sigma_c), ncol = trueC_len)

      Yt = true_val + Xt %*% C + ut
      Yc = Xc %*% C + uc

    }
  }else if(d == 2){

    if(X_dist == "equi"){

      n2 <- as.integer(sqrt(n))^2  # sqrt(n) is not an integer
      n3 <- n - n2 # to make the sample size back to n

      grid_pt <- seq(from = -1, to = 1, length.out = sqrt(n2))
      X_2 <- data.matrix(expand.grid(grid_pt, grid_pt))
      X <- rbind(X_2, matrix(0, nrow = n3, ncol = 2))

      tind <- X[,1] < 0 & X[,2] < 0
    }

    Xt <- X[tind == 1, ,drop = F]
    Xc <- X[tind == 0, ,drop = F]
    sigma_t <- sig[tind == 1]
    sigma_c <- sig[tind == 0]

    if(spec == "Lin"){

      ut <- rep(stats::rnorm(length(sigma_t), mean = 0, sd = sigma_t), trueC_len)
      ut <- matrix(ut, nrow = length(sigma_t), ncol = trueC_len)
      uc <- rep(stats::rnorm(length(sigma_c), mean = 0, sd = sigma_c), trueC_len)
      uc <- matrix(uc, nrow = length(sigma_c), ncol = trueC_len)

      Yt = true_val + Xt[, 1, drop = F] %*% C + Xt[, 2, drop = F] %*% C + ut
      Yc = Xc[, 1, drop = F] %*% C + Xc[, 2, drop = F] %*% C + uc

    }
  }

  res <- list(Yt = Yt, Yc = Yc, Xt = Xt, Xc = Xc, sigma_t = sigma_t, sigma_c = sigma_c)
  return(res)
}

#' CI Generation Corresponding to the Specified Method
#'
#' @param met CI generation method to use,
#' among \code{c("Ex", "Csvtv", "Ex_mm", "Csvtv_mm", "rdr")} so far.
#' @inheritParams c_hat_lower_RD
#' @param alpha a desired level of non-coverage.
#' @param rdr_swap \code{TRUE} if swapped value of \code{rdrobust} result will be used;
#' set \code{TRUE} if \code{x < 0} is treated.
#' @param trueC a vector of true smoothness parameters used in simulations.
#' @param C_len the number of smoothness parameters to adapt to.
#' @param Csvtv_const the constant which governs the degree of conservativeness of the
#' adaptive and minimax procedures.
#' @param lower \code{TRUE} in order to generate a lower one-sided confidence interval
#' @param rdr_met which row of the \code{rdrobust} result matrix will be used; the default is 3.
#'
#' @return a vector of lower and upper ends of the CI.
#' @export
#'
#' @examples obs_data <- rdadapt::gen_obs(100, 1, 1/2, "Lin", "het", "equi")
#' Xt <- obs_data$Xt
#' Xc <- obs_data$Xc
#' sigma_t <- obs_data$sigma_t
#' sigma_c <- obs_data$sigma_c
#' Yt <- obs_data$Yt
#' Yc <- obs_data$Yc
#' trueC <- (1:5)/5
#' C_len <- 5
#' CI_gen_met("Ex", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05, trueC, C_len)
#' CI_gen_met("Ex_2", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05, trueC, C_len)
#' CI_gen_met("rdr", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05)
#' CI_gen_met("Ex_mm", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05, trueC, C_len,
#' Csvtv_const = 1, lower = TRUE)
#' CI_gen_met("Nomon", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05, trueC, C_len)
#' CI_gen_met("Csvtv_hbrd", Xt, Xc, c(1), sigma_t, sigma_c, Yt, Yc, 0.05, trueC, C_len,
#' Csvtv_const = 2, lower = FALSE)
CI_gen_met<- function(met = c("Ex", "Csvtv", "Ex_2", "Csvtv_2", "Ex_mm", "Csvtv_mm",
                              "Csvtv_hbrd", "Nomon", "rdr"),
                      Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                      trueC = NULL, C_len = NULL, Csvtv_const = NULL, lower = FALSE,
                      rdr_swap = TRUE, rdr_met = 3){

  met = match.arg(met)

  if(met == "Ex"){

    C_min <- min(trueC)
    C_max <- max(trueC)
    Cvec <- seq(from = C_min, to = C_max, length.out = C_len)

    CI <- CI_adpt(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, lower)

  }else if(met == "Csvtv"){

    C_min <- min(trueC)
    C_max <- Csvtv_const * max(trueC)
    Cvec <- seq(from = C_min, to = C_max, length.out = C_len)

    CI <- CI_adpt(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, lower)

  }else if(met == "Ex_2"){

    C_min <- min(trueC)
    C_max <- max(trueC)
    Cvec <- seq(from = C_min, to = C_max, length.out = 2)

    CI <- CI_adpt(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, lower)

  }else if(met == "Csvtv_2"){

    C_min <- min(trueC)
    C_max <- Csvtv_const * max(trueC)
    Cvec <- seq(from = C_min, to = C_max, length.out = 2)

    CI <- CI_adpt(Cvec, Cbar, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha, lower)

  }else if(met == "Ex_mm"){

    C_max <- max(trueC)

    if(lower == FALSE){

      CI <- CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, sigma_t, sigma_c, alpha)

    }else{

      CI_l <- c_hat_lower_RD(stats::qnorm(1 - alpha), C_max, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c,
                         Yt, Yc, alpha)
      CI <- c(CI_l, Inf)
    }


  }else if(met == "Nomon"){

    C_max <- max(trueC)

    if(lower == FALSE){

      CI <- CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, NULL, sigma_t, sigma_c, alpha)

    }else{

      CI_l <- c_hat_lower_RD(stats::qnorm(1 - alpha), C_max, C_max, Xt, Xc, NULL, sigma_t, sigma_c,
                             Yt, Yc, alpha)
      CI <- c(CI_l, Inf)
    }


  }else if(met == "Csvtv_mm"){

    C_max <- Csvtv_const * max(trueC)

    if(lower == FALSE){

      CI <- CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, sigma_t, sigma_c, alpha)

    }else{

      CI_l <- c_hat_lower_RD(stats::qnorm(1 - alpha), C_max, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c,
                             Yt, Yc, alpha)
      CI <- c(CI_l, Inf)
    }

  }else if(met == "Csvtv_hbrd"){

    C_min <- min(trueC)
    C_max <- Csvtv_const * max(trueC)
    Cvec <- seq(from = C_min, to = C_max, length.out = C_len)

    CI_1 <- CI_adpt(Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha/2, lower)

    if(lower == FALSE){

      CI_2 <- CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, sigma_t, sigma_c, alpha/2)

    }else{

      CI_l <- c_hat_lower_RD(stats::qnorm(1 - alpha/2), C_max, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c,
                             Yt, Yc, alpha/2)
      CI_2 <- c(CI_l, Inf)
    }

    CI <- c(max(CI_1[1], CI_2[1]), min(CI_1[2], CI_2[2]))

  }else if(met == "rdr"){

    X <- as.numeric(rbind(Xt, Xc))
    Y <- c(Yt, Yc)


    if(lower == FALSE){

      rdr_res_all <- rdrobust::rdrobust(y = Y, x = X, level = (1 - alpha) * 100)

      rdr_res <- rdr_res_all$ci
      CI <- c(rdr_res[3, 1], rdr_res[3, 2])
      if(rdr_swap) CI <- -CI[c(2:1)]

    }else{

      rdr_res_all <- rdrobust::rdrobust(y = Y, x = X, level = (1 - 2 * alpha) * 100)

      rdr_res <- rdr_res_all$ci
      CI_all <- c(rdr_res[3, 1], rdr_res[3, 2])
      CI <- c(-Inf, CI_all[2])
      if(rdr_swap){

        CI_all <- -CI_all[c(2:1)]
        CI <- c(CI_all[1], Inf)
      }
    }
  }
  return(CI)
}

#' CI Result Generation
#'
#' @param n sample size.
#' @param spec_list a list generated by \code{spec_set}.
#' @param method_sym a vector CI construction methods to be used; currently supports
#' \code{c("Ex", "rdr")}.
#' @param true_val the true parameter value; default is 1. It can be a vector for
#' power function calculations.
#' @param power generates power function result if \code{TRUE}.
#' @inheritParams CI_gen_met
#'
#' @return a vector of (coverage, CI_length) for each method - a vector of
#' length \code{len(method_sym) * 2}.
#' @export
#'
#' @examples n <- 100
#' spec_list <- rdadapt::spec_set(1, 2, 0.5)
#' method_sym <- c("Ex", "rdr")
#' res_gen(n, spec_list, method_sym, 5, 2)
#' res_gen(n, spec_list, method_sym, 5, 2, FALSE, c(0, 1), TRUE)
res_gen <- function(n, spec_list, method_sym, C_len, Csvtv_const, lower = FALSE, true_val = 1,
                    power = FALSE){

  if(power == FALSE){

    d <- spec_list$d
    mon_ind <- spec_list$mon_ind
    alpha <- spec_list$alpha
    trueC <- spec_list$trueC
    trueC_min <- min(trueC)
    trueC_max <- max(trueC)
    fcn_spec <- spec_list$fcn_spec
    sd_spec <- spec_list$sd_spec
    X_dist <- spec_list$X_dist

    obs_data <- gen_obs(n, d, trueC, fcn_spec, sd_spec, X_dist , true_val)
    Xt <- obs_data$Xt
    Xc <- obs_data$Xc

    sigma_t <- obs_data$sigma_t
    sigma_c <- obs_data$sigma_c
    Yt <- obs_data$Yt
    Yc <- obs_data$Yc

    m_len <- length(method_sym)
    trueC_len <- length(trueC)

    res <- numeric(2 * trueC_len * m_len)

    for(j in 1:trueC_len){

      Yt_j <- as.numeric(Yt[, j])
      Yc_j <- as.numeric(Yc[, j])

      res_j <- numeric(2 * m_len)

      for(i in 1:m_len){

        met <- method_sym[i]
        ind <- 2 * (i - 1) + 1

        CI <- CI_gen_met(met, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt_j, Yc_j, alpha,
                         trueC, C_len, Csvtv_const, lower)

        res_j[ind] <- as.numeric(CI[1] < true_val & CI[2] > true_val)

        if(lower == FALSE){
          res_j[ind + 1] <- CI[2] - CI[1]
        }else{
          res_j[ind + 1] <- true_val - CI[1]
        }

      }

      ind_1 <- 2 * m_len * (j - 1) + 1
      ind_2 <- 2 * m_len * j

      res[ind_1:ind_2] <- res_j

    }
  }else{

    d <- spec_list$d
    mon_ind <- spec_list$mon_ind
    alpha <- spec_list$alpha
    trueC <- spec_list$trueC
    trueC_min <- min(trueC)
    trueC_max <- max(trueC)
    fcn_spec <- spec_list$fcn_spec
    sd_spec <- spec_list$sd_spec
    X_dist <- spec_list$X_dist
    trueC_sc <- spec_list$trueC_sc

    true_val_len <- length(true_val)
    m_len <- length(method_sym)

    res <- numeric(true_val_len * m_len)

    for(j in 1:true_val_len){

      true_val_j <- true_val[j]
      obs_data <- gen_obs(n, d, trueC_sc, fcn_spec, sd_spec, X_dist, true_val_j)
      Xt <- obs_data$Xt
      Xc <- obs_data$Xc

      sigma_t <- obs_data$sigma_t
      sigma_c <- obs_data$sigma_c
      Yt <- obs_data$Yt
      Yc <- obs_data$Yc

      res_j <- numeric(m_len)

      for(i in 1:m_len){

        met <- method_sym[i]

        CI <- CI_gen_met(met, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, alpha,
                         trueC, C_len, Csvtv_const, lower)

        res_j[i] <- as.numeric(CI[1] > 0)

      }

      ind_1 <- m_len * (j - 1) + 1
      ind_2 <- m_len * j

      res[ind_1:ind_2] <- res_j

    }
  }


  return(res)

}


#' Dataframe Generation After Simulations
#'
#' @param res_all data matrix generated from simulations.
#' @param trueC a vector of true smoothness parameters used.
#' @param method_name a vector of names of CI construction methods used.
#'
#' @return a list containing coverage and length results stored as \code{data.frame}.
#' @export
#'
#' @examples trueC <- c(0.1, 0.2)
#' method_name <- c("a", "b")
#' m_len <- length(method_name)
#' trueC_len <- length(trueC)
#' nIters <- 50
#' rowlen <- m_len * trueC_len * 2
#' res_all <- matrix(stats::rnorm(rowlen * nIters), nrow = rowlen, ncol = nIters)
#' res_form(res_all, trueC, method_name)
res_form <- function(res_all, trueC, method_name){

  m_len <- length(method_name)
  trueC_len <- length(trueC)

  resMean <- rowMeans(res_all)
  col_C <- rep(trueC, each = m_len)
  col_met <- rep(method_name, trueC_len)
  col_cov <- resMean[c(TRUE, FALSE)]
  col_len <- resMean[c(FALSE, TRUE)]

  cov_data <- data.frame(C = col_C, cov_prob = col_cov, method = col_met)
  len_data <- data.frame(C = col_C, avg_len = col_len, method = col_met)

  res <- list(cov_data = cov_data, len_data = len_data)

  return(res)
}


#' Dataframe Generation After Simulations (Power Function)
#'
#' @param true_val_vec a vector of true values.
#' @inheritParams res_form
#'
#' @return a list containing power results.
#' @export
#'
#' @examples true_val_vec <- c(0.1, 0.2)
#' method_name <- c("a", "b")
#' m_len <- length(method_name)
#' true_val_len <- length(true_val_vec)
#' nIters <- 50
#' rowlen <- m_len * true_val_len
#' res_all <- matrix(stats::rnorm(rowlen * nIters), nrow = rowlen, ncol = nIters)
#' res_form(res_all, true_val_vec, method_name)
res_form_pow <- function(res_all, true_val_vec, method_name){

  m_len <- length(method_name)
  true_val_len <- length(true_val_vec)

  resMean <- rowMeans(res_all)
  col_tv <- rep(true_val_vec, each = m_len)
  col_met <- rep(method_name, true_val_len)

  pow_data <- data.frame(true_val = col_tv, rej_prob = resMean, method = col_met)

  return(pow_data)
}
