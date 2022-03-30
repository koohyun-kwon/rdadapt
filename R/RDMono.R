#' Inference in RD under Monotonicity
#'
#' Description to be supplied.
#'
#' Under development.
#'
#' @param X a vector or a matrix of the running variables, with \code{ncol(X)} being
#' the dimension of the turnning variable when \code{X} is a matrix.
#' @param Y a vector of the outcome variable.
#' @param t.ind  a vector of treatment indicators, with \code{t.ind = 1} indicating
#' treated observations and \code{t.ind = 0} control observations.
#' @param C bound on the first derivative of the regression function; if \code{C} is missing,
#' \code{C = Inf} is used when \code{method} corresponds to one of
#' \code{c("adpt.L", "adpt.U")} (see below).
#' @param mon.ind index for monotone variables, a subset of \code{1:ncol(X)} if \code{X}
#' is a matrix, either 1 or 0 if \code{X} is a vector.
#' @param C.l lower bound for the Lipschitz constants to which the adaptive one-sided CI
#' adapt to.
#' @param C.u upper bound for the Lipschitz constants to which the adaptive one-sided CI
#' adapt to.
#' @param method type of confidence interval to be used. The options are:
#'
#' \describe{
#'
#' \item{"mm"}{Minimax two-sided confidence interval}
#'
#' \item{"adpt.L"}{Adaptive one-sided lower confidence interval}
#'
#' \item{"adpt.U"}{Adaptive one-sided upper confidence interval}
#'
#' }
#' @param se.method method for estimating the standard error of the estimate. (to be edited later)
#' @param se.init method for estimating the initial variance for computing the optimal bandwidth.
#' (to be edited later)
#' @param alpha determines confidence level, \code{1 - alpha}.
#' @param N number of nearest neighbors to be used when "nn" is specified in \code{se.method}.
#'
#' @return a list of following components
#' \describe{
#'
#' \item{\code{ci}}{Confidence interval}
#'
#' \item{\code{ht, hc}}{Bandwidths used for treated and control observations}
#'
#' }
#' @export
#'
#' @examples print("under development")
RDMono <- function(X, Y, t.ind, C, mon.ind, C.l, C.u, method = c("mm", "adpt.L", "adpt.U"),
                   se.method, se.init, alpha, N){


}
