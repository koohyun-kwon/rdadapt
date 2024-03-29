% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_fcn_RD.R
\name{Norm}
\alias{Norm}
\title{(Weighted) \eqn{L_p} Norm}
\usage{
Norm(X, p = 1, invw = 1)
}
\arguments{
\item{X}{A numeric matrix.}

\item{p}{Order of the norm; the default is \code{p = 1} (L1 norm).}

\item{invw}{Inverse weights for each component,
with the length equal to either 1 or \code{ncol(X)}; the default is
\code{invw = 1}.}
}
\value{
A vector of the dimension equal to the \code{nrow(X)}.
}
\description{
Calculates rowwise \eqn{||x||_p} for each row \eqn{x}
of a matrix \code{X}.
}
\details{
\code{X} should be a matrix, not a vector.
}
\examples{
X <- matrix(c(1, -2, -3, 4, 5, -6), nrow = 3, ncol = 2)
Norm(X)
weights <- c(1, 2)
Norm(X = X, invw = weights)
}
