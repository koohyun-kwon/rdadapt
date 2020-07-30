#' Minimax CI length for regression function at a point
#'
#' @inheritParams invmod_RD
#' @param C a scalar smoothness parameter.
#' @param alpha the desired level of non-coverage probability.
#'
#' @return the length of CI corresponding to the modulus value \code{b}
#' @export
#'
#' @examples n <- 500
#' d <- 2
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] > 0 & X[, 2] > 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' C <- 1/2
#' mon_ind <- c(1, 2)
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' CI_length_RD(1, C, Xt, Xc, mon_ind, sigma_t , sigma_c, 0.05)
CI_length_RD <- function(b, C, Xt, Xc, mon_ind, sigma_t , sigma_c, alpha) {

  om_inv <- invmod_RD(b, rep(C, 2), Xt, Xc, mon_ind, sigma_t, sigma_c)

  bc <- om_inv$bc
  bt <- om_inv$bt
  om_inv_t <- om_inv$delta_t
  om_inv_c <- om_inv$delta_c
  om_inv <- sqrt(om_inv_t^2 + om_inv_c^2)

  if (om_inv == 0) return(Inf)

  Kt <- K_fun(bt, rep(C, 2), Xt, mon_ind)
  gf_ip_iota_t <- sum(bt * Kt / sigma_t^2)

  Kc <- K_fun(bc, rep(C, 2), Xc, mon_ind)
  gf_ip_iota_c <- sum(bc * Kc / sigma_c^2)

  sd <- sqrt((om_inv_t/gf_ip_iota_t)^2 + (om_inv_c/gf_ip_iota_c)^2)

  bias <- .5 * (b - (om_inv^2 /  (gf_ip_iota_t + gf_ip_iota_c)))
  cva <- ifelse(abs(bias / sd) > 3,
                abs(bias / sd) + stats::qnorm(1 - alpha),
                sqrt(stats::qchisq(1 - alpha, df = 1, ncp = (bias / sd)^2)))

  return(2 * cva * sd)
}

#' Minimax Confidence Interval
#'
#' @inheritParams c_hat_lower_RD
#' @param C_max the worst-case smoothness parameter.
#' @param alpha the desired level of non-coverage
#' @param opt_b provided if the optimal modulus value is known; default is \code{NULL}.
#' @param min_half_length provided if the optimal half-length is known;
#' default is \code{NULL}.
#' @param maxb.const governs the optimization range; default is 10.
#' @param Prov.Plot if \code{TRUE}, provides a plot that can be used to check the optimization
#' worked well; default is \code{FALSE}.
#'
#' @return the left and right ends of the confidence interval
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
#' CI_minimax_RD(Yt, Yc, Xt, Xc, C_max, mon_ind, sigma_t, sigma_c, 0.05)
CI_minimax_RD <- function(Yt, Yc, Xt, Xc, C_max, mon_ind, sigma_t, sigma_c,
                          alpha, opt_b = NULL, min_half_length = NULL,
                          maxb.const = 10, Prov.Plot = FALSE) {

  if(is.null(opt_b) | is.null(min_half_length)){

    # modres <- modsol_RD(0,rep(C_max,2), Xt, Xc, mon_ind,
    #                     sigma_t, sigma_c)
    # minbt <- modres$bt
    # minbc <- modres$bc
    # minb <- minbt + minbc

    minbt <- minb_fun(rep(C_max, 2), Xt, mon_ind)
    minbc <- minb_fun(rep(C_max, 2), Xc, mon_ind, swap = T)
    minb <- minbt + minbc

    modres_2 <- modsol_RD(stats::qnorm(1 - alpha/2), rep(C_max,2), Xt, Xc,
                          mon_ind, sigma_t, sigma_c)

    maxb <- maxb.const * (modres_2$bt + modres_2$bc)

    CI_length_sol <- stats::optimize(CI_length_RD, interval = c(minb, maxb),
                                     C = C_max, Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                                     sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)

    min_half_length <- CI_length_sol$objective / 2
    opt_b <- CI_length_sol$minimum

    if(Prov.Plot == TRUE){

      numgrid = 100
      xintv = seq(from = minb, to = maxb, length.out = numgrid)
      yvec = numeric(numgrid)
      for(i in 1:numgrid){

        yvec[i] = CI_length_RD(xintv[i], C = C_max, Xt = Xt, Xc = Xc, mon_ind = mon_ind,
                               sigma_t = sigma_t, sigma_c = sigma_c, alpha = alpha)
      }

      graphics::plot(xintv, yvec, type = "l", xlab = "modulus", ylab = "CI_length")
      graphics::abline(v = opt_b, col = "red", lty = 2)
    }

  }

  invmod_opt <- invmod_RD(opt_b, rep(C_max, 2), Xt, Xc, mon_ind, sigma_t, sigma_c)
  ht <- invmod_opt$bt
  hc <- invmod_opt$bc

  opt_Lhat <- Lhat_fun_RD(0, C_max, C_max, Xt, Xc, mon_ind,
                          sigma_t, sigma_c, Yt, Yc, ht, hc) +
              a_fun(0, C_max, C_max, Xt, Xc, mon_ind, sigma_t, sigma_c, ht, hc)

  res <- c(opt_Lhat - min_half_length, opt_Lhat + min_half_length)

  return(res)
}
