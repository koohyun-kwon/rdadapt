#' High-level CI Generation
#'
#' Generates a confidence interval according to a given method.
#'
#' This is used for MC simulations.
#'
#' @param CI.method method name used to generate the confidence interval.
#' @param xt running variable values for the treated.
#' @param xc running variable values for the control.
#' @param yt outcome variable values for the treated.
#' @param yc outcome variable values for the control.
#' @param sig.t standard deviations for the treated.
#' @param sig.c standard deviations for the control.
#' @param spec.etc list of other specification parameters.
#' @param C.one a single Lipschitz coefficient to adapt to;
#'  needed for \code{CI.method = "adpt.one"}.
#'
#' @return lower and upper ends of a CI, returned as a list
#' (\code{ci.l} and \code{ci.u}).
#' @export
#'
CI.gen <- function(CI.method = c("mm.smallC", "mm.largeC", "RDH", "adpt", "RDR",
                                 "RDR.L", "adpt.one"),
                   xt, xc, yt, yc, sig.t, sig.c, spec.etc, C.one){

  CI.method = match.arg(CI.method)

  alpha <- spec.etc$alpha
  mon_ind <- spec.etc$mon_ind
  se.method <- spec.etc$se.method
  se.init <- spec.etc$se.init
  t.dir <- spec.etc$t.dir

  if(CI.method == "adpt"){

    C_l <- spec.etc$C.l
    C_u <- spec.etc$C.u
    C <- spec.etc$C

    res.adpt <- CI_adpt_opt(C_l, C_u, C, xt, xc, mon_ind, yt, yc, alpha, se.method,
                            sig.t, sig.c, se.init = se.init, t.dir = t.dir)
    res <- list(ci = res.adpt, sd = -1) # "sd" is not well defined for an adaptive procedure

  }else if(CI.method == "mm.smallC" | CI.method == "mm.largeC"){

    if(CI.method == "mm.smallC"){

      C_max <- spec.etc$C.small
    }else{

      C_max <- spec.etc$C.large
    }

    res.mm <- CI_minimax_RD(yt, yc, xt, xc, C_max, mon_ind, se.method, se.init,
                            t.dir, alpha, sigma_t = sig.t, sigma_c = sig.c)

    res <- list(ci = res.mm$ci, sd = res.mm$sd)

  }else if(CI.method == "RDH"){

    M.RDH <- spec.etc$M.RDH
    se.initial <- spec.etc$se.initial.RDH
    se.method <- spec.etc$se.method.RDH

    x <- c(as.vector(xt), as.vector(xc)) # for the case where xt and xc are 1-column matrices
    y <- c(yt, yc)
    sig2 <- c(sig.t, sig.c)^2

    if(se.initial == "supplied.var"){

      # 1st stage: supplied; 2nd stage: se.method
      d <- RDHonest::RDData(data.frame(y, x, sigma2 = sig2), 0)
      res.RDH <- RDHonest::NPRHonest.fit(d, M = M.RDH, opt.criterion = "FLCI",
                                         alpha = alpha, sclass = "H", se.method = se.method)

    }else{

      # 1st stage: se.initial; 2nd stage: nn
      d <- RDHonest::RDData(data.frame(y, x), 0)
      res.RDH <- RDHonest::NPRHonest.fit(d, M = M.RDH, opt.criterion = "FLCI",
                                         alpha = alpha, sclass = "H", se.initial = se.initial)
    }


    ci.res <- c(res.RDH$estimate - res.RDH$hl, res.RDH$estimate + res.RDH$hl)
    ci.res <- -c(ci.res[2], ci.res[1]) # RDHonest treats x > 0 as treated

    res <- list(ci = ci.res, sd = res.RDH$sd)
  }else if(CI.method == "RDR"){

    x <- c(as.vector(xt), as.vector(xc)) # for the case where xt and xc are 1-column matrices
    y <- c(yt, yc)

    rdr.res <- rdrobust::rdrobust(y, x, level = 100 * (1 - alpha))
    # the third row contains the robust CI result
    ci.res <- -c(rdr.res$ci[3, 2], rdr.res$ci[3, 1]) # rdrobust treats x > 0 as treated
    res <- list(ci = ci.res, sd = rdr.res$se[3, 1])

  }else if(CI.method == "RDR.L"){

    x <- c(as.vector(xt), as.vector(xc)) # for the case where xt and xc are 1-column matrices
    y <- c(yt, yc)

    rdr.res <- rdrobust::rdrobust(y, x)
    # rdrobust treats x > 0 as treated
    ci.res <- -c(rdr.res$coef[3, 1] + rdr.res$se[3, 1] * stats::qnorm(1 - alpha), -Inf)
    res <- list(ci = ci.res, sd = rdr.res$se[3, 1])

  }else if(CI.method == "adpt.one"){

    C_l <- C.one
    C_u <- C.one
    C <- spec.etc$C

    res.adpt <- CI_adpt_opt(C_l, C_u, C, xt, xc, mon_ind, yt, yc, alpha, se.method,
                            sig.t, sig.c, se.init = se.init, t.dir = t.dir)
    res <- list(ci = res.adpt, sd = -1)
  }

  return(res)
}
