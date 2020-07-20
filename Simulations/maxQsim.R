#' Simulation to confirm whether the quntile of "LHS" does not
#' depend on the value of tau

# Design
n <- 500
d <- 2
X <- matrix(rnorm(n * d), nrow = n, ncol = d)
tind <- X[, 1] > 0 & X[, 2] > 0
Xt <- X[tind == 1, ,drop = F]
Xc <- X[tind == 0, ,drop = F]
Cvec <- seq(from = 0.2, to = 1, by = 0.2)
mon_ind <- c(1, 2)
sigma <- rnorm(n)^2 + 1
sigma_t <- sigma[tind == 1]
sigma_c <- sigma[tind == 0]

# Simulation

taulen <- 50
tauvec <- seq(from = 0.05 / length(Cvec), to = 0.3, length.out = taulen)
resvec1 <- numeric(taulen)
resvec2 <- numeric(taulen)

for(i in 1:taulen){

  delta <- qnorm(1 - tauvec[i])
  covmat <- cov_mat_calc(delta, Cvec, Xt, Xc, mon_ind, sigma_t, sigma_c)
  resvec1[i] <- max_Q(covmat, 0.05, num_sim = 10^4)
  resvec2[i] <- max_Q(covmat, 0.05, num_sim = 10^5)
  print(i)
}

plot(tauvec, resvec, type = "l")
title(main = "n = 500, d = 2, num_sim = (10^4, 10^5)")
lines(tauvec, resvec2, col = "red")
