#' Controls the Iterations in the bfsl Algorithm
#'
#' \code{bfsl_control} allows the user to set some characteristics of the \code{bfsl}
#' best-fit straight line algorithm.
#'
#' @param tol A positive numeric value specifying the tolerance level for the
#' convergence criterion
#' @param maxit A positive integer specifying the maximum number of iterations allowed.
#'
#' @return A \code{list} with two components named as the arguments.
#'
#' @seealso \code{\link{bfsl}}
#'
#' @examples bfsl_control(tol = 1e-8, maxit = 1000)
#'
#' @export
bfsl_control = function(tol = 1e-10, maxit = 100) {
  list(tol = tol, maxit = maxit)
}

#' Calculates the Best-fit Straight Line
#'
#' \code{bfsl} calculates the best-fit straight line to independent points with
#' (possibly correlated) normally distributed errors in both coordinates.
#'
#' \code{bfsl} provides the general least-squares estimation solution to the
#' problem of fitting a straight line to independent data with (possibly
#' correlated) normally distributed errors in both \code{x} and \code{y}.
#'
#' With \code{sd_x = 0} the (weighted) ordinary least squares solution is
#' obtained. The calculated standard errors of the slope and intercept
#' multiplied with \code{sqrt(chisq)} correspond to the ordinary least squares
#' standard errors.
#'
#' With \code{sd_x = c}, \code{sd_y = d}, where \code{c} and \code{d} are
#' positive numbers, and \code{r = 0} the Deming regression solution is obtained.
#' If additionally \code{c = d}, the orthogonal distance regression solution,
#' also known as major axis regression, is obtained.
#'
#' Setting \code{sd_x = sd(x)}, \code{sd_y = sd(y)} and \code{r = 0} leads to
#' the geometric mean regression solution, also known as reduced major
#' axis regression or standardised major axis regression.
#'
#' The goodness of fit metric \code{chisq} is a weighted reduced chi-squared
#' statistic. It compares the deviations of the points from the fit line to the
#' assigned measurement error standard deviations. If \code{x} and \code{y} are
#' indeed related by a straight line, and if the assigned measurement errors
#' are correct (and normally distributed), then \code{chisq} will equal 1. A
#' \code{chisq > 1} indicates underfitting: the fit does not fully capture the
#' data or the measurement errors have been underestimated. A \code{chisq < 1}
#' indicates overfitting: either the model is improperly fitting noise, or the
#' measurement errors have been overestimated.
#'
#' @param x A vector of \emph{x} observations or a data frame (or an
#' object coercible by \code{\link{as.data.frame}} to a data frame) containing
#' the named vectors \emph{x}, \emph{y}, and optionally \emph{sd_x},
#' \emph{sd_y} and \emph{r}. If weights \emph{w_x} and \emph{w_y} are given,
#' then \emph{sd_x} and \emph{sd_y} are calculated from \emph{sd_x = 1/sqrt(w_x)}
#' and \emph{sd_y = 1/sqrt(w_y)}. Specifying \code{y}, \code{sd_x}, \code{sd_y}
#' or \code{r} directly as function arguments overwrites these variables in the
#' data structure.
#' @param y A vector of \emph{y} observations.
#' @param sd_x A vector of \emph{x} measurement error standard
#' deviations. If it is of length one, all data points are assumed to have the
#' same \emph{x} standard deviation.
#' @param sd_y A vector of \emph{y} measurement error standard
#' deviations. If it is of length one, all data points are assumed to have the
#' same \emph{y} standard deviation.
#' @param r A vector of correlation coefficients between errors in
#' \emph{x} and \emph{y}. If it is of length one, all data points are assumed to
#' have the same correlation coefficient.
#' @param control A list of control settings. See \code{\link{bfsl_control}}
#' for the names of the settable control values and their effect.
#'
#' @return An object of class "\code{bfsl}", which is a \code{list} containing
#' the following components:
#' \item{coefficients}{A \code{2x2} matrix with columns of the fitted coefficients
#' (intercept and slope) and their standard errors.}
#' \item{chisq}{The goodness of fit  (see Details).}
#' \item{control}{The control \code{list} used, see the \code{control} argument.}
#' \item{convInfo}{A \code{list} with convergence information.}
#' \item{call}{The matched call.}
#' \item{data}{A \code{list} containing \code{x}, \code{y}, \code{sd_x}, \code{sd_y}
#' and \code{r}.}
#'
#' @references York, D. (1968). Least squares fitting of a straight line with
#' correlated errors. \emph{Earth and Planetary Science Letters}, 5, 320–324,
#' https://doi.org/10.1016/S0012-821X(68)80059-7
#'
#' @examples
#' x = pearson_york$x
#' y = pearson_york$y
#' sd_x = 1/sqrt(pearson_york$w_x)
#' sd_y = 1/sqrt(pearson_york$w_y)
#' bfsl(x, y, sd_x, sd_y)
#'
#' fit = bfsl(pearson_york)
#' plot(fit)
#'
#' @export
bfsl = function(x, y = NULL, sd_x = 0, sd_y = 1, r = 0, control = bfsl_control()) {

  # dispatch variables if first argument is a data frame, list, etc.
  if (!is.vector(x) || is.list(x)) {
    df = as.data.frame(x)
    if (exists("x", df)) {
      x = df$x
    } else {
      stop("x variable is missing in data structure.")
    }
    if (missing("y") && exists("y", df)) {
      y = df$y
    }
    if (missing("sd_x") && exists("sd_x", df)) {
      sd_x = df$sd_x
    }
    if (missing("sd_y") && exists("sd_y", df)) {
      sd_y = df$sd_y
    }
    if (missing("sd_x") && exists("w_x", df)) {
      sd_x = 1/sqrt(df$w_x)
    }
    if (missing("sd_y") && exists("w_y", df)) {
      sd_y = 1/sqrt(df$w_y)
    }
    if (missing("r") && exists("r", df)) {
      r = df$r
    }
  }

  # record the function call
  cl = match.call()

  # iterations control
  ctrl = bfsl_control()
  if(!missing(control)) {
    control = as.list(control)
    ctrl[names(control)] = control
  }

  # check arguments
  if (is.null(y)) { stop("Argument 'y' is missing.") }
  if (!is.numeric(x) || !is.numeric(y) || !is.vector(x) || !is.vector(y) ||
      length(x) != length(y)) {
    stop("Arguments 'x' and 'y' must be numeric vectors of equal length.")
  }
  n = length(x)
  if (!is.numeric(sd_x) || !is.vector(sd_x) || (length(sd_x) != n &&
                                                length(sd_x) != 1)) {
    stop("Argument 'sd_x' must be a numeric vector the same length as 'x' or of length 1.")
  }
  if (!is.numeric(sd_y) || !is.vector(sd_y) || (length(sd_y) != n &&
                                                length(sd_y) != 1)) {
    stop("Argument 'sd_y' must be a numeric vector the same length as 'y' or of length 1.")
  }
  if (!is.numeric(r) || !is.vector(r) || (length(r) != n && length(r) != 1)) {
    stop("Argument 'r' must be a numeric vector the same length as 'x' or of length 1.")
  }

  # ordinary least squares for initial guess value of slope
  x0 = as.matrix(cbind(1,x))
  b = (solve(t(x0) %*% x0) %*% t(x0) %*% y)[2]

  # expand sd_x and sd_y if they are of length 1
  if (length(sd_x) == 1) {
    sd_x = rep(sd_x, n)
  }
  if (length(sd_y) == 1) {
    sd_y = rep(sd_y, n)
  }

  # weights
  sd_x = pmax(sd_x, .Machine$double.eps)
  sd_y = pmax(sd_y, .Machine$double.eps)
  wX = 1/sd_x^2
  wY = 1/sd_y^2
  alpha = sqrt(wX*wY)

  # iterate until the slope converges to within the tolerance
  b_diff = 2*control$tol
  i = 0
  while (b_diff > control$tol && i < control$maxit) {
    i =  i + 1
    b0 = b
    W = wX*wY/(wX + b^2*wY - 2*b*r*alpha)
    meanX = sum(W*x)/sum(W)
    meanY = sum(W*y)/sum(W)
    U = x - meanX
    V = y - meanY
    beta = W*(U/wY + b*V/wX - (b*U + V)*(r/alpha))
    b = sum(W*beta*V)/sum(W*beta*U)  # slope
    b_diff = abs(b - b0)
  }

  a = meanY - b*meanX  # intercept
  meanx = sum(W*(meanX + beta))/sum(W)
  u = (meanX + beta) - meanx
  sd_b = sqrt(1/(sum(W*u^2)))  # standard error of b
  sd_a = sqrt(1/sum(W) + meanx^2*sd_b^2)  # standard error of a
  chisq = sum(W*(y - b*x - a)^2)/(length(x)-2)  # goodness of fit

  # convergence information
  convInfo = list(isConv = (i<control$maxit), finIter = i, finTol = b_diff)

  # fitted coefficients including standard errors
  coef = c(a, b)
  names(coef) = c("Intercept", "Slope")
  sterr = c(sd_a, sd_b)
  coefficients = cbind(coef, sterr)
  dimnames(coefficients) = list(names(coef), c("Estimate", "Std. Error"))

  bfsl.out = list(coefficients = coefficients, chisq = chisq,
                  control = control, convInfo = convInfo, call = cl,
                  data = list(x = x, y = y, sd_x = sd_x, sd_y = sd_y, r = r))

  class(bfsl.out) = "bfsl"

  return(bfsl.out)
}


#' Print Method for bfsl Results
#'
#' \code{print} method for class "\code{bfsl}".
#'
#' @param x An object of class "\code{bfsl}".
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to \code{print.default}.
#'
#' @export
print.bfsl = function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Best-fit straight line\n\n")

  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE, ...)

  cat("\nGoodness of fit:\n")
  cat(format(x$chisq, digits = digits))

  invisible(x)
}

#' Plot Method for bfsl Results
#'
#' Plot method for objects of class "\code{bfsl}".
#'
#' \code{plot.bfsl} plots the data points with error bars and the calculated
#' best-fit straight line.
#'
#' @param x An object of class "\code{bfsl}".
#' @param grid If \code{TRUE} (default) grid lines are plotted.
#' @param ... Further parameters to be passed to the plotting routines.
#'
#' @importFrom graphics abline plot points segments
#'
#' @export
plot.bfsl = function(x, grid = TRUE, ...)
{
  obj = x
  x = obj$data$x
  y = obj$data$y
  sd_x = obj$data$sd_x
  sd_y = obj$data$sd_y

  bw_x = 0.01*diff(range(x))  # bar width in x direction
  bw_y = 0.015*diff(range(y))  # bar width in y direction

  plot(x, y, xlim = c(min(x-sd_x), max(x+sd_x)),
       ylim = c(min(y-sd_y), max(y+sd_y)), type = "n", ...)
  if (grid) grid()
  points(x, y, ...)

  # error bars
  segments(x, y-sd_y, x, y+sd_y)
  segments(x-bw_x, y+sd_y, x+bw_x, y+sd_y)
  segments(x-bw_x, y-sd_y,x+bw_x, y-sd_y)
  segments(x-sd_x, y,x+sd_x, y)
  segments(x+sd_x, y-bw_y, x+sd_x, y+bw_y)
  segments(x-sd_x, y-bw_y, x-sd_x, y+bw_y)

  # fit line
  abline(coef = obj$coefficients[,1])
}
