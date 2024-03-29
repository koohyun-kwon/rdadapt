% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_intsec_RD.R
\name{CI_adpt}
\alias{CI_adpt}
\title{Adaptive Confidence Interval}
\usage{
CI_adpt(
  Cvec,
  Cbar,
  Xt,
  Xc,
  mon_ind,
  sigma_t,
  sigma_c,
  Yt,
  Yc,
  alpha,
  lower = TRUE,
  num_sim = 10^5,
  delta_init = 1.96,
  bmat_init_L,
  bmat_L,
  tau_res
)
}
\arguments{
\item{Cvec}{a sequence of smoothness parameters}

\item{Cbar}{the Lipschitz coefficient for the largest function space we consider}

\item{Xt}{\eqn{n_t} by \eqn{k} design matrix for the treated units.}

\item{Xc}{\eqn{n_c} by \eqn{k} design matrix for the control units.}

\item{mon_ind}{index number for monotone variables.}

\item{sigma_t}{standard deviation of the error term for the treated units
(either length 1 or \eqn{n_t}).}

\item{sigma_c}{standard deviation of the error term for the control units
(either length 1 or \eqn{n_c}).}

\item{Yt}{outcome value for the treated group observations.}

\item{Yc}{outcome value for the control group observations.}

\item{alpha}{desired upper quantile value.}

\item{lower}{calculate a lower one-sided confidence interval if \code{TRUE};
calculate a two-sided CI otherwise.}

\item{num_sim}{number of simulations used to calculate the quantile;
the default is \code{10^5}.}

\item{delta_init}{the value of \eqn{\delta} to be used in simulating the quantile;
theoretically, its value does not matter asymptotically. Its default value is 1.96.}

\item{bmat_init_L}{the matrix of modulus values corresponding to
\code{delta_init} and \code{Cvec}; it can be left unspecified.}

\item{bmat_L}{the matrix of modulus values corresponding to
the optimal \eqn{\delta} and \code{Cvec}; it can be left unspecified.}

\item{tau_res}{outcome value from the function \code{tau_res}; it can be
left unspecified.}
}
\value{
confidence interval endpoints.
}
\description{
Calculates the one-sided or two-sided adaptive confidence interval
}
\details{
This returns either one-sided lower CI or two-sided Bonferroni CI using
the function CI_adpt_L; it seems this function will be substituted by
\code{CI_adpt_opt} defined below.
}
\examples{
n <- 500
d <- 2
X <- matrix(rnorm(n * d), nrow = n, ncol = d)
tind <- X[, 1] < 0 & X[, 2] < 0
Xt <- X[tind == 1, ,drop = FALSE]
Xc <- X[tind == 0, ,drop = FALSE]
mon_ind <- c(1, 2)
sigma <- rnorm(n)^2 + 1
sigma_t <- sigma[tind == 1]
sigma_c <- sigma[tind == 0]
Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
CI_adpt((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = FALSE)
CI_adpt((1:5)/5, 1, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = TRUE)
CI_adpt((1:5)/5, Inf, Xt, Xc, mon_ind, sigma_t, sigma_c, Yt, Yc, 0.05, lower = TRUE)
}
