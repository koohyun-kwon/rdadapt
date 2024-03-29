% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCres_comb.R
\name{f.eval}
\alias{f.eval}
\title{Regression Function Evaluation}
\usage{
f.eval(which.fun, x, C, th)
}
\arguments{
\item{which.fun}{a positive integer; curretly supports 1:6.}

\item{x}{vector of regressor values.}

\item{C}{smoothness parameter.}

\item{th}{true RD parameter.}
}
\value{
vector of evaluated values with length \code{length(x)}
}
\description{
Calculates the regression function value given a specification.
}
\details{
When \code{which.fun} is 5 or 6, \code{x} should be a matrix with two columns.
}
\examples{
f.eval(2, runif(100,-1,1), 1, 1)
}
