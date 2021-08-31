#' Function \eqn{f_1}
#'
#' Function \eqn{f_1} in the version 20201104
#'
#' @param x vector of regressor values.
#' @param C smoothness parameter.
#' @param th true RD parameter.
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x <- runif(100, -1, 1)
#' f.1(x, 1, 1)
f.1 <- function(x, C, th){

  C * x + th * (x < 0)

}


#' Function \eqn{f_2}
#'
#' Function \eqn{f_2} in the version 20201104
#'
#' @inheritParams f.1
#' @param b1 the first knot value; the default is \code{b1 = 1/3}.
#' @param b2 the second knot value; the default is \code{b2 = 2/3}.
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x <- runif(100, -1, 1)
#' f.2(x, 1, 1)
f.2 <- function(x, C, th, b1 = 1/3, b2 = 2/3){

  f.2c <- function(x, C){

    x.new <- x * (x >= 0)
    1.5 * C * (x.new^2 - 2 * (x.new - b1)^2 * (x.new > b1) +
                 2 * (x.new - b2)^2 * (x.new > b2))
  }

  (x < 0) * (-f.2c(-x, C) + th) + (x >= 0) * f.2c(x, C)
}

#' Function \eqn{f_3}
#'
#' Function \eqn{f_3} in the version 20201104
#'
#' @inheritParams f.1
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x <- runif(100, -1, 1)
#' f.3(x, 1, 1)
f.3 <- function(x, C, th){

  (C / 4) * (x^3 + x) + th *(x < 0)
}

#' Function \eqn{f_4}
#'
#' Function \eqn{f_4} in the version 20201105
#'
#' @inheritParams f.1
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x <- runif(100, -1, 1)
#' f.4(x, 1, 1)
f.4 <- function(x, C, th){

  f.4c <- function(x, C){

    x.new <- x * (x >= 0)
    C * ((3 * x.new + 1)^(1/3) - 1)
  }

  (x < 0) * (-f.4c(-x, C) + th) + (x >= 0) * f.4c(x, C)
}


#' Function \eqn{f_5}
#'
#' Function \eqn{f_5} in the version 20201104
#'
#' @param x.mat n by 2 matrix
#' @inheritParams f.1
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x.mat <- cbind(runif(100, -1, 1), runif(100, -1, 1))
#' f.5(x.mat, 1, 1)
f.5 <- function(x.mat, C, th){

  C * (x.mat[, 1] + x.mat[, 2]) + th * (x.mat[, 1] < 0 & x.mat[, 2] < 0)
}

#' Function \eqn{f_6}
#'
#' Function \eqn{f_6} in the version 20201104
#'
#' @inheritParams f.5
#'
#' @return vector of regression function values
#' @export
#'
#' @examples x.mat <- cbind(runif(100, -1, 1), runif(100, -1, 1))
#' f.6(x.mat, 1, 1)
f.6 <- function(x.mat, C, th){

  f.6c <- function(x.mat, C){

    x1.new <- x.mat[, 1] * (x.mat[, 1] >= 0 | x.mat[, 2] >= 0)
    x2.new <- x.mat[, 2] * (x.mat[, 1] >= 0 | x.mat[, 2] >= 0)

    C * (sqrt(2 * (x1.new + x2.new) + 3) - sqrt(3))
  }

  (x.mat[, 1] < 0 & x.mat[, 2] < 0) * (-f.6c(-x.mat, C) + th) +
    (x.mat[, 1] >= 0 | x.mat[, 2] >= 0) * f.6c(x.mat, C)
}
