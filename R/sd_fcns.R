#' Nearest Neighborhood Estimation of Conditional Standard Deviation
#'
#' Standard deviation is calculated by the formula (23) of Armstrong and Kolesár (2020)
#'
#' Currently, this is not being used, although this works for general dimensions;
#' debugging to be done.
#'
#' @inheritParams CI_adpt
#' @param N the number of nearest neighbors; default is \code{N = 3}.
#'
#' @return a list containing conditional standard deviation estimates for treated
#' observations (\code{sigma.t}) and control observations (\code{sigma.c})
#'
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and honest confidence intervals in nonparametric regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' }
#'
#' @export
#'
sigmaNN2 <- function(Xt, Xc, Yt, Yc, N = 3){ # N = 3 might not be suitable for d > 1 case

  # Euclidean distance between each observations based on X's
  # Returns a n_q by n_q matrix for q = {t,c}
  dist.t <- as.matrix(stats::dist(Xt, diag = TRUE, upper = TRUE))
  dist.c <- as.matrix(stats::dist(Xc, diag = TRUE, upper = TRUE))

  # Prevent choosing itself
  diag(dist.t) <- Inf
  diag(dist.c) <- Inf

  nt <- nrow(Xt)
  nc <- nrow(Xc)
  yt.nb <- matrix(0, nt, N)
  yc.nb <- matrix(0, nc, N)

  for(i in 1:nt){

    order.i <- order(dist.t[i, ])
    yt.nb[i, ] <- Yt[order.i[1:N]]
  }

  for(i in 1:nc){

    order.i <- order(dist.c[i, ])
    yc.nb[i, ] <- Yc[order.i[1:N]]
  }

  diff.t <- Yt - rowMeans(yt.nb)
  diff.c <- Yc - rowMeans(yc.nb)

  # Formula (23) of Armstrong and Kolesár (2020)
  sigma.t <- sqrt((N / (N + 1)) * diff.t^2)
  sigma.c <- sqrt((N / (N + 1)) * diff.c^2)

  return(list(sigma.t = sigma.t, sigma.c = sigma.c))
}


#' Nearest Neighborhood Estimation of Conditional Standard Deviation
#'
#' Standard deviation is calculated by the formula (23) of Armstrong and Kolesár (2020)
#'
#' For now, this function works only for one-dimensional cases.
#'
#' @inheritParams CI_adpt
#' @param N the number of nearest neighbors; default is \code{N = 3}.
#' @param t.dir treatment direction; \code{t.dir = "left"} if \eqn{x < 0} is treated.
#' Otherwise, \code{t.dir = "right"}.
#'
#' @return a list containing conditional standard deviation estimates for treated
#' observations (\code{sigma.t}) and control observations (\code{sigma.c})
#'
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and honest confidence intervals in nonparametric regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' }
#'
#' @export
#'
#' @examples n <- 500
#' d <- 1
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' sigma <- rnorm(n)^2 + 1
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' sigmaNN(Xt, Xc, Yt, Yc, t.dir = "left")
sigmaNN <- function(Xt, Xc, Yt, Yc, t.dir = c("left", "right"), N = 3){

  t.dir <- match.arg(t.dir)

  RDH.sigmaNN <- function (X, Y, J, weights = rep(1L, length(X))){

    n <- length(X)
    sigma2 <- matrix(nrow = n, ncol = NCOL(Y)^2)
    for (k in seq_along(X)) {
      s <- max(k - J, 1):k
      d <- sort(abs(c(X[s[-length(s)]], X[k:min(k + J, n)][-1]) -
                      X[k]))[J]
      ind <- (abs(X - X[k]) <= d)
      ind[k] <- FALSE
      Jk <- sum(weights[ind])
      sigma2[k, ] <- Jk/(Jk + weights[k]) * if (NCOL(Y) >
                                                1)
        as.vector(outer(Y[k, ] - colSums(weights[ind] *
                                           Y[ind, ])/Jk, Y[k, ] - colSums(weights[ind] *
                                                                            Y[ind, ])/Jk))
      else (Y[k] - sum(weights[ind] * Y[ind])/Jk)^2
    }
    drop(sigma2)
  }

  y <- c(Yt, Yc)
  x <- c(as.vector(Xt), as.vector(Xc))
  d <- RDHonest::RDData(data.frame(y = y, x = x), 0)

  if(t.dir == "left"){
    # Flip p and m since RDHonest assumes x > 0 is treated
    sigma.t <- RDH.sigmaNN(d$Xm, d$Ym, N)
    sigma.c <- RDH.sigmaNN(d$Xp, d$Yp, N)
  }else{

    sigma.t <- RDH.sigmaNN(d$Xp, d$Yp, N)
    sigma.c <- RDH.sigmaNN(d$Xm, d$Ym, N)
  }

  return(list(sigma.t = sqrt(sigma.t), sigma.c = sqrt(sigma.c)))
}

#' Standard Deviation Calculation by Silverman's Rule
#'
#' Standard deviation is calculated using the first-stage bandwidth chosen by
#' Silverman's Rule, as in RDHonest::NPRPrelimVar.fit.
#'
#' This works only for cases with one-dimensional running variables.
#'
#' @inheritParams CI_adpt
#' @param t.dir treatment direction; \code{t.dir = "left"} if
#' \eqn{x < 0} is treated. Otherwise, \code{t.dir = "right"}.
#'
#' @return a list containing conditional standard deviation estimates for treated
#' observations (\code{sigma.t}) and control observations (\code{sigma.c})
#' @export
#'
#' @examples X <- matrix(rnorm(500 * 1), nrow = 500, ncol = 1)
#' tind <- X[, 1] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' sigma <- rep(1, 500)
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' sigmaSvm(Xt, Xc, Yt, Yc, "left")
sigmaSvm <- function(Xt, Xc, Yt, Yc, t.dir = c("left", "right")){

  Xt <- as.vector(Xt)
  Xc <- as.vector(Xc)
  X <- c(Xt, Xc)
  Y <- c(Yt, Yc)

  # Minimum value of meaningful bandwidth
  hmin <- max(sort(abs(unique(Xt)))[2], sort(abs(unique(Xc)))[2])

  h1 <- max(1.84 * stats::sd(X)/sum(length(X))^(1/5), hmin)

  # First-stage nonparametric regression
  # order = 0, since we don't assume bounded 2nd derivative
  drf <- data.frame(y = Y, x = X)
  d <- RDHonest::RDData(drf, 0)
  r1 <- RDHonest::NPRreg.fit(d = d, h = h1, order = 0, se.method = "EHW")

  r1$sigma2p <- r1$sigma2p * length(r1$sigma2p)/(length(r1$sigma2p) - 1)
  r1$sigma2m <- r1$sigma2m * length(r1$sigma2m)/(length(r1$sigma2m) - 1)

  if(t.dir == "left"){

    # Flip p and m since RDHonest assumes x > 0 is treated
    sigma.c <- rep(mean(r1$sigma2p), length(d$Xp))
    sigma.t <- rep(mean(r1$sigma2m), length(d$Xm))
  }else{

    sigma.t <- rep(mean(r1$sigma2p), length(d$Xp))
    sigma.c <- rep(mean(r1$sigma2m), length(d$Xm))
  }

  return(list(sigma.t = sqrt(sigma.t), sigma.c = sqrt(sigma.c)))
}


#' Standard Deviation Calculation by Silverman's Rule (Tesing Version for Multi-dim)
#'
#' Standard deviation is calculated using Silverman's Rule,
#' as in RDHonest::NPRPrelimVar.fit
#'
#' This works only for one-dimensional cases.
#'
#' @inheritParams CI_adpt
#'
#' @return a list containing conditional standard deviation estimates for treated
#' observations (\code{sigma.t}) and control observations (\code{sigma.c})
#' @export
#'
#' @examples n <- 500
#' d <- 1
#' X <- matrix(rnorm(n * d), nrow = n, ncol = d)
#' tind <- X[, 1] < 0
#' Xt <- X[tind == 1, ,drop = FALSE]
#' Xc <- X[tind == 0, ,drop = FALSE]
#' sigma <- rep(1, n)
#' sigma_t <- sigma[tind == 1]
#' sigma_c <- sigma[tind == 0]
#' Yt = 1 + rnorm(length(sigma_t), mean = 0, sd = sigma_t)
#' Yc = rnorm(length(sigma_c), mean = 0, sd = sigma_c)
#' sigmaSvm.test(Xt, Xc, Yt, Yc)
sigmaSvm.test <- function(Xt, Xc, Yt, Yc){

  # Convert running variable data into matrices (for dim = 1)
  if(!is.matrix(Xt)) Xt <- matrix(Xt, ncol = 1)
  if(!is.matrix(Xc)) Xc <- matrix(Xc, ncol = 1)

  d <- ncol(Xt)
  n <- nrow(Xt) + nrow(Xc)

  # Minimum bandwidth configuration
  Xmin.vec <- numeric(d)
  for(i in 1:d){
    Xt.i <- Xt[, i]
    Xc.i <- Xc[, i]

    Xmin.vec[i] <- max(sort(abs(unique(Xt.i)))[2], sort(abs(unique(Xc.i)))[2])
  }

  X <- rbind(Xt, Xc)
  Y <- c(Yt, Yc)

  # Calculating the initial bandwidth using multivariate kernel density rule of thumb
  # Based on Table 6.3 and (6.43) of Scott (1992)
  # When d = 1, equivalent to h1 defined in Section 4.2 of Imbens and Kalyanaraman (2012),
  # except for the minimum bandwidth configuration

  h.vec <- numeric(d)
  for(i in 1:d){

    X.i <- X[, i]
    h.vec[i] <- max(1.74 * (4 / ((d + 2) * n))^(1/(d + 4)) * stats::sd(X.i), Xmin.vec[i])
  }

  # Local constant regression
  if(d == 1){
    bwres.t <- np::npregbw(Yt ~ as.vector(Xt), bws = h.vec, bandwidth.compute = FALSE)
    npreg.res.t <- np::npreg(Yt ~ as.vector(Xt), bws = bwres.t, ckertype = "uniform",
                           newdata = data.frame(Xt = 0))
    bwres.c <- np::npregbw(Yc ~ as.vector(Xc), bws = h.vec, bandwidth.compute = FALSE)
    npreg.res.c <- np::npreg(Yc ~ as.vector(Xc), bws = bwres.c, ckertype = "uniform",
                             newdata = data.frame(Xc = 0))

    resid.t <- (Yt - npreg.res.t$mean)[abs(Xt[, 1]) < h.vec[1]]
    resid.c <- (Yc - npreg.res.c$mean)[abs(Xc[, 1]) < h.vec[1]]
  }else if(d == 2){
    Xt1 <- Xt[, 1]
    Xt2 <- Xt[, 2]
    Xc1 <- Xc[, 1]
    Xc2 <- Xc[, 2]
    bwres.t <- np::npregbw(Yt ~ Xt1 + Xt2, bws = h.vec, bandwidth.compute = FALSE)
    npreg.res.t <- np::npreg(Yt ~ Xt1 + Xt2, bws = bwres.t, ckertype = "uniform",
                             newdata = data.frame(Xt1 = 0, Xt2 = 0))
    bwres.c <- np::npregbw(Yc ~ Xc1 + Xc2, bws = h.vec, bandwidth.compute = FALSE)
    npreg.res.c <- np::npreg(Yc ~ Xc1 + Xc2, bws = bwres.c, ckertype = "uniform",
                             newdata = data.frame(Xc1 = 0, Xc2 = 0))

    resid.t <- (Yt - npreg.res.t$mean)[abs(Xt[, 1]) < h.vec[1] & abs(Xt[, 2]) < h.vec[2]]
    resid.c <- (Yc - npreg.res.c$mean)[abs(Xc[, 1]) < h.vec[1] & abs(Xc[, 2]) < h.vec[2]]
  }

  # Estimation of conditional variance under homoskedasticity

  sigma2t <- resid.t^2 * length(resid.t) / (length(resid.t) - d)
  sigma2c <- resid.c^2 * length(resid.c) / (length(resid.c) - d)

  sigma.t <- rep(sqrt(mean(sigma2t)), nrow(Xt))
  sigma.c <- rep(sqrt(mean(sigma2c)), nrow(Xc))

  return(list(sigma.t = sigma.t, sigma.c = sigma.c))
}
