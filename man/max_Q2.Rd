% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_intsec_RD.R
\name{max_Q2}
\alias{max_Q2}
\title{Quantile of Maximum of Multivariate Normal Random Variable 2}
\usage{
max_Q2(covmat, alpha)
}
\arguments{
\item{covmat}{covariance matrix of the multivariate normal random variable.}

\item{alpha}{desired upper quantile value.}
}
\value{
\eqn{(1 - \alpha)}th quantile of maximum of a multivariate normal random variable.
}
\description{
Calculates the quantile of maximum of a multivariate normal random variable,
using multivariate normal distribution function calculation.
}
\details{
This function was created for a test purpose; it is not used since it is slow for
large \eqn{J}.
}
\examples{
covmat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
max_Q2(covmat, 0.05)
}
