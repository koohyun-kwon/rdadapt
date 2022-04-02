#' Quantile of Folded Normal Distribution
#'
#' Calculates the \eqn{1 - \alpha} quantile of \eqn{|N(t, 1)|}.
#'
#' @param t non-centrality parameter.
#' @param alpha probability
#'
#' @return the quantile value
#' @export
#'
#' @examples cv_a(1, 0.05)
cv_a <- function(t, alpha){

  q_sq <- stats::qchisq(1 - alpha, 1, ncp = t^2)
  res <- sqrt(q_sq)

  return(res)
}

#' Minimum Half-length Function
#'
#' Calculates the minimum half-length, as a function of the modulus value.
#'
#' In the paper, the minimum half-length is expressed as a function of
#' \eqn{\delta};
#' here, we use the relationship that \eqn{\delta} is the inverse modulus value
#' given \eqn{b}.
#' The returned value will be minimized with respect to \eqn{b} later
#'
#' @inheritParams invmod_RD
#' @param C a scalar smoothness parameter.
#' @param alpha the desired level of non-coverage probability.
#'
#' @return the minimum half-length corresponding to the modulus value \code{b}
#' @export
#'
#' @examples X <- matrix(rnorm(500 * 2), nrow = 500, ncol = 2)
#' tind <- X[, 1] > 0 & X[, 2] > 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' sigma <- rnorm(500)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' CI_length_RD(1.1, 1/2, Xt, Xc, c(1, 2), sigma_t , sigma_c, 0.05)
CI_length_RD <- function(b, C, Xt, Xc, mon_ind, sigma_t , sigma_c, alpha) {

  # Calculate delta corresponding to b
  # This is because sup_bias calculation is with respect to delta
  om_inv <- invmod_RD(b, rep(C, 2), Xt, Xc, mon_ind, sigma_t, sigma_c)
  om_inv_t <- om_inv$delta_t
  om_inv_c <- om_inv$delta_c
  delta <- sqrt(om_inv_t^2 + om_inv_c^2)

  # Retrieve b_c and b_t to speed up the computation of sup_bias
  bc <- om_inv$bc
  bt <- om_inv$bt

  # Deal with the case where delta = 0; this can happen when b is very small
  # Basically, this is the case where bandwidth is too small (no effective obs)
  # We prevent this case by letting the function return Inf
  if ((om_inv_t + om_inv_c) == 0) return(Inf)

  # Remaining codes are for the case with delta > 0
  sup_bias <- sup_bias_Lhat_RD(delta, C, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                               bt, bc)
  sup_bias <- sup_bias - a_fun(delta, C, C, Xt, Xc, mon_ind, sigma_t, sigma_c,
                               bt, bc)
  sd <- sd_Lhat_RD(delta, C, C, Xt, Xc, mon_ind, sigma_t, sigma_c, bt, bc)

  t <- sup_bias / sd
  res <- cv_a(t, alpha) * sd

  return(res)
}

#' Minimax Confidence Interval
#'
#' Calculates the minimax confidence interval.
#'
#' So far, conditional variance estimation works only for one-dimensional case.
#'
#' @inheritParams c_hat_lower_RD
#' @param se.method the standard deviation estimation method.
#' @param se.init the standard deviation estimation method for choosing an optimal estimator.
#' @param sigma_t supplied variance for treated observations.
#' @param sigma_c supplied variance for control observations.
#' @param sigma_t.init supplied first-stage variance for treated observations.
#' @param sigma_c.init supplied first-stage variance for control observations.
#' @param C_max the worst-case smoothness parameter.
#' @param t.dir treatment direction; \code{t.dir = "left"} if \eqn{x < 0} is treated.
#' Otherwise, \code{t.dir = "right"}. This should specified only for one-dimensional cases.
#' @param alpha the desired level of non-coverage
#' @param N the number of neighbors to be used when \code{se.method = "nn"};
#' the default is \code{N = 3}.
#' @param opt_b provided if the optimal modulus value is known; default is \code{NULL}.
#' @param min_half_length provided if the optimal half-length is known;
#' default is \code{NULL}.
#' @param maxb.const governs the optimization range; default is 10.
#' @param Prov.Plot if \code{TRUE}, provides a plot that can be used to check the optimization
#' worked well; default is \code{FALSE}.
#' @param len.return if \code{TRUE}, returns only the optimal half-length;
#' default is \code{FALSE}.
#'
#' @return returns a list with the confidence interval (\code{ci}), the standard deviation
#' of the estimator (\code{sd}), and the bandwidths used for the treated observations and
#' the control observations (\code{h.t} and \code{h.c})
#'
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
#' Yt <- 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc <- rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' C_max <- 1
#' CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, "nn.test", "S.test", alpha = 0.05)
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
#' CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, "nn", "Silverman", t.dir = "left",
#' alpha = 0.05)
#' CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, "nn", "nn", t.dir = "left",
#' alpha = 0.05)
CI_minimax_RD <- function(Yt, Yc, Xt, Xc, C_max, mon_ind,
                          se.method = c("nn", "supplied", "nn.test"),
                          se.init = c("Silverman", "nn", "supplied", "supp.sep",
                                      "S.test"),
                          t.dir = c("left", "right"),
                          alpha, N = 3, sigma_t, sigma_c,
                          sigma_t.init, sigma_c.init,
                          opt_b = NULL, min_half_length = NULL, maxb.const = 10,
                          Prov.Plot = FALSE, len.return = FALSE) {

  # For case with d = 1, Xt & Xc might be vector, not matrix
  if(!is.matrix(Xt)) Xt <- matrix(Xt, ncol = 1)
  if(!is.matrix(Xc)) Xc <- matrix(Xc, ncol = 1)

  se.method = match.arg(se.method)
  se.init = match.arg(se.init)
  t.dir = match.arg(t.dir)

  # The sorting process is necessary due to the way variance estimation
  # functions work (Only supports d = 1 case)
  if(ncol(Xt) == 1){
    sorted.t <- sort(Xt, index.return = T)
    sorted.c <- sort(Xc, index.return = T)
    Xt <- matrix(sorted.t$x, ncol = 1)
    Xc <- matrix(sorted.c$x, ncol = 1)

    # Y indexes should match the ones for X
    Yt <- Yt[sorted.t$ix]
    Yc <- Yc[sorted.c$ix]

    # Sigma indexes should match the ones for X
    if(!missing(sigma_t) & !missing(sigma_c)){

      sigma_t <- sigma_t[sorted.t$ix]
      sigma_c <- sigma_c[sorted.c$ix]
    }
    if(!missing(sigma_t.init) & !missing(sigma_c.init)){

      sigma_t.init <- sigma_t.init[sorted.t$ix]
      sigma_c.init <- sigma_c.init[sorted.c$ix]
    }
  }

  if(se.init != "supplied"){
    # Check whether this part is documented somewhere
    if(se.init == "Silverman"){

      if(ncol(Xt) > 1) stop("Multi-dimension not supported for now")

      sigma.init <- sigmaSvm(Xt, Xc, Yt, Yc, t.dir)
      sigma_t.init <- sigma.init$sigma.t
      sigma_c.init <- sigma.init$sigma.c

    }else if(se.init == "nn"){

      if(ncol(Xt) > 1) stop("Multi-dimension not supported for now")

      # Currently testing whether it works (so is the test done?)
      sigma.init <- sigmaNN(Xt, Xc, Yt, Yc, t.dir)
      sigma_t.init <- sigma.init$sigma.t
      sigma_c.init <- sigma.init$sigma.c

    # What is this doing? - Silverman's rule for multiple variable
    # The variance estimation function is still under test
    }else if(se.init == "S.test"){

      sigma.init <- sigmaSvm.test(Xt, Xc, Yt, Yc)
      sigma_t.init <- sigma.init$sigma.t
      sigma_c.init <- sigma.init$sigma.c
    }
  }

  if(is.null(opt_b) | is.null(min_half_length)){

    # Minimum modulus calculation for determining the optimization range
    minbt <- minb_fun(rep(C_max, 2), Xt, mon_ind)
    minbc <- minb_fun(rep(C_max, 2), Xc, mon_ind, swap = T)
    minb <- minbt + minbc

    # Determining the upper end of the optimization range by heuristics
    modres_U <- modsol_RD(stats::qnorm(1 - alpha/2), rep(C_max,2), Xt, Xc,
                          mon_ind, sigma_t.init, sigma_c.init)
    maxb <- maxb.const * (modres_U$bt + modres_U$bc)

    CI_length_sol <- stats::optimize(CI_length_RD, interval = c(minb, maxb),
                                     C = C_max, Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                                     sigma_t = sigma_t.init, sigma_c = sigma_c.init,
                                     alpha = alpha)
    min_half_length.init <- CI_length_sol$objective
    opt_b <- CI_length_sol$minimum

    # This provides a way to check if the optimization was done correctly
    if(Prov.Plot == TRUE){

      numgrid = 100
      xintv = seq(from = minb, to = maxb, length.out = numgrid)
      yvec = numeric(numgrid)
      for(i in 1:numgrid){

        yvec[i] = CI_length_RD(xintv[i], C = C_max, Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                               sigma_t = sigma_t.init, sigma_c = sigma_c.init,
                               alpha = alpha)
      }

      graphics::plot(xintv, yvec, type = "l", xlab = "modulus", ylab = "CI_length")
      graphics::abline(v = opt_b, col = "red", lty = 2)
    }
  }

  if(len.return == TRUE){

    return(min_half_length.init) # returns the value without the second-stage sd estimation

  }else{

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

    # Retrieving the optimal value of delta; the first stage sd values are used
    invmod_opt <- invmod_RD(opt_b, rep(C_max, 2), Xt, Xc, mon_ind,
                            sigma_t.init, sigma_c.init)
    ht <- invmod_opt$bt
    hc <- invmod_opt$bc
    del_t <- invmod_opt$delta_t
    del_c <- invmod_opt$delta_c
    delta <- sqrt(del_t^2 + del_c^2)

    # Calculates the form of the optimal estimator; the first stage sd values are used
    optL.res <- Lhat_fun_RD(delta, C_max, C_max, Xt, Xc, mon_ind,
                            sigma_t.init, sigma_c.init, Yt, Yc, ht, hc,
                            ret.w = TRUE)
    opt_Lhat <- optL.res$est -
      a_fun(delta, C_max, C_max, Xt, Xc, mon_ind, sigma_t.init, sigma_c.init,
            ht, hc)

    # Re-calculates the half-length using the second stage sd estimation
    sd <- sd_w(optL.res$w_t, optL.res$w_c, sigma_t.new, sigma_c.new)

    min_half_length <- sd * min_half_length.init /
      sd_Lhat_RD(delta, C_max, C_max, Xt, Xc, mon_ind, sigma_t.init,
                 sigma_c.init, ht, hc)

    # Returns a list
    res <- list(ci = c(opt_Lhat - min_half_length, opt_Lhat + min_half_length),
                sd = sd, h.t = ht/C_max, h.c = hc/C_max)

    return(res)
  }
}
