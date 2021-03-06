% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_mm.R
\name{CI_minimax_RD_len}
\alias{CI_minimax_RD_len}
\title{Length of Minimax Confidence Interval}
\usage{
CI_minimax_RD_len(
  Xt,
  Xc,
  C_max,
  mon_ind,
  sigma_t,
  sigma_c,
  alpha,
  opt_b = NULL,
  min_half_length = NULL,
  maxb.const = 10,
  Prov.Plot = FALSE
)
}
\arguments{
\item{Xt}{\eqn{n_t} by \eqn{k} design matrix for the treated units.}

\item{Xc}{\eqn{n_c} by \eqn{k} design matrix for the control units.}

\item{C_max}{the worst-case smoothness parameter.}

\item{mon_ind}{index number for monotone variables.}

\item{sigma_t}{standard deviation of the error term for the treated units
(either length 1 or \eqn{n_t}).}

\item{sigma_c}{standard deviation of the error term for the control units
(either length 1 or \eqn{n_c}).}

\item{alpha}{the desired level of non-coverage}

\item{opt_b}{provided if the optimal modulus value is known; default is \code{NULL}.}

\item{min_half_length}{provided if the optimal half-length is known;
default is \code{NULL}.}

\item{maxb.const}{governs the optimization range; default is 10.}

\item{Prov.Plot}{if \code{TRUE}, provides a plot that can be used to check the optimization
worked well; default is \code{FALSE}.}
}
\value{
the minimax CI half-length value.
}
\description{
Same as \code{CI_minimax_RD}, except that this function only returns the minimax length.
}
\examples{
n <- 500
d <- 2
X <- matrix(rnorm(n * d), nrow = n, ncol = d)
tind <- X[, 1] > 0 & X[, 2] > 0
Xt <- X[tind == 1, ,drop = FALSE]
Xc <- X[tind == 0, ,drop = FALSE]
mon_ind <- c(1, 2)
sigma <- rnorm(n)^2 + 1
sigma_t <- sigma[tind == 1]
sigma_c <- sigma[tind == 0]
C_max <- 1
CI_minimax_RD_len(Xt, Xc, C_max, mon_ind, sigma_t, sigma_c, 0.05)
}
