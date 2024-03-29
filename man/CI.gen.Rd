% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_gen_all.R
\name{CI.gen}
\alias{CI.gen}
\title{High-level CI Generation}
\usage{
CI.gen(
  CI.method = c("mm.smallC", "mm.largeC", "RDH", "adpt", "RDR", "RDR.L", "adpt.one"),
  xt,
  xc,
  yt,
  yc,
  sig.t,
  sig.c,
  spec.etc,
  C.one
)
}
\arguments{
\item{CI.method}{method name used to generate the confidence interval.}

\item{xt}{running variable values for the treated.}

\item{xc}{running variable values for the control.}

\item{yt}{outcome variable values for the treated.}

\item{yc}{outcome variable values for the control.}

\item{sig.t}{standard deviations for the treated.}

\item{sig.c}{standard deviations for the control.}

\item{spec.etc}{list of other specification parameters.}

\item{C.one}{a single Lipschitz coefficient to adapt to;
needed for \code{CI.method = "adpt.one"}.}
}
\value{
lower and upper ends of a CI, returned as a list
(\code{ci.l} and \code{ci.u}).
}
\description{
Generates a confidence interval according to a given method.
}
\details{
This is used for MC simulations.
}
