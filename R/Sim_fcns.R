#' Sample Generation for Simulations
#'
#' @param n sample size.
#' @param d number of running variables.
#' @param C_pair a pair of smoothness parameters; \code{C_pair[1] <= C_pair[2]} should hold.
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
#' @examples gen_obs(100, 1, c(1/2, 1), "Const", "het", "equi")
gen_obs <- function(n, d, C_pair, spec = c("Const","Lin"), sd_spec = c("hom", "het"),
                    X_dist = c("equi"), true_val = 1){

  X_dist <- match.arg(X_dist)
  sd_spec <- match.arg(sd_spec)
  spec <- match.arg(spec)

  if(sd_spec == "hom"){

    sig = rep(1, n)

  }else if(sd_spec == "het"){

    sig <- stats::rnorm(n, 0, 1/sqrt(2))^2 + 1/2
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

    if(spec == "Const"){

      Yt = true_val + stats::rnorm(length(sigma_t), mean = 0, sd = sigma_t)
      Yc = stats::rnorm(length(sigma_c), mean = 0, sd = sigma_c)

    }else if(spec == "Lin"){

      C_min <- C_pair[1]
      Yt = true_val + C_min * Xt + stats::rnorm(length(sigma_t), mean = 0, sd = sigma_t)
      Yc = C_min * Xc + stats::rnorm(length(sigma_c), mean = 0, sd = sigma_c)
    }
  }

  res <- list(Yt = Yt, Yc = Yc, Xt = Xt, Xc = Xc, sigma_t = sigma_t, sigma_c = sigma_c)
  return(res)
}
