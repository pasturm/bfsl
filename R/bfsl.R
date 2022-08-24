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
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class "\code{bfsl}", which is a \code{list} containing
#' the following components:
#' \item{coefficients}{A \code{2x2} matrix with columns of the fitted coefficients
#' (intercept and slope) and their standard errors.}
#' \item{chisq}{The goodness of fit  (see Details).}
#' \item{fitted.values}{The fitted mean values.}
#' \item{residuals}{The residuals, that is \code{y} observations minus fitted values.}
#' \item{df.residual}{The residual degrees of freedom.}
#' \item{cov.ab}{The covariance of the slope and intercept.}
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
#' x = pearson_york_data$x
#' y = pearson_york_data$y
#' sd_x = 1/sqrt(pearson_york_data$w_x)
#' sd_y = 1/sqrt(pearson_york_data$w_y)
#' bfsl(x, y, sd_x, sd_y)
#' bfsl(y~x, pearson_york_data, sd_x, sd_y)
#'
#' fit = bfsl(pearson_york_data)
#' plot(fit)
#'
#' @export
bfsl <- function(...) { UseMethod("bfsl") }

#' @rdname bfsl
#' @export
bfsl.default = function(x, y = NULL, sd_x = 0, sd_y = 1, r = 0,
                        control = bfsl_control(), ...) {

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

  out = bfsl_fit(x, y, sd_x, sd_y, r, ctrl, cl)

  return(out)
}

#' @param formula A formula specifying the bivariate model (as in \code{\link{lm}},
#' but here only \code{y ~ x} makes sense).
#' @param data A data.frame containing the variables of the model.
#'
#' @rdname bfsl
#' @export
bfsl.formula = function(formula, data = parent.frame(), sd_x, sd_y, r = 0,
                        control = bfsl_control(), ...) {

  if (missing("sd_x") && exists("sd_x", data)) {
    sd_x = data$sd_x
  }
  if (missing("sd_y") && exists("sd_y", data)) {
    sd_y = data$sd_y
  }
  if (missing("sd_x") && exists("w_x", data)) {
    sd_x = 1/sqrt(data$w_x)
  }
  if (missing("sd_y") && exists("w_y", data)) {
    sd_y = 1/sqrt(data$w_y)
  }
  if (missing("r") && exists("r", data)) {
    r = data$r
  }

  cl = match.call()
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels = TRUE
  mf[[1L]] = quote(stats::model.frame)
  mf = eval(mf, parent.frame())

  y = mf[,1]
  x = mf[,2]

  # iterations control
  ctrl = bfsl_control()
  if(!missing(control)) {
    control = as.list(control)
    ctrl[names(control)] = control
  }

  out = bfsl_fit(x, y, sd_x, sd_y, r, ctrl, cl)

  return(out)
}

bfsl_fit = function(x, y, sd_x, sd_y, r, control, cl) {

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

  fitted.values = b*x+a
  residuals = y - fitted.values
  df.residual = length(y)-2
  cov.ab = -meanX*sd_b^2

  bfsl.out = list(coefficients = coefficients, chisq = chisq,
                  fitted.values = fitted.values, residuals = residuals,
                  df.residual = df.residual, cov.ab = cov.ab,
                  control = control, convInfo = convInfo, call = cl,
                  data = list(x = x, y = y, sd_x = sd_x, sd_y = sd_y, r = r))

  class(bfsl.out) = "bfsl"

  return(bfsl.out)
}

#' Print Method for bfsl Results
#'
#' \code{print} method for class \code{"bfsl"}.
#'
#' @param x An object of class "\code{bfsl}".
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to \code{print.default}.
#'
#' @export
print.bfsl = function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE, ...)

  cat("\n")
  invisible(x)
}

#' Plot Method for bfsl Results
#'
#' \code{plot.bfsl} plots the data points with error bars and the calculated
#' best-fit straight line.
#'
#' @param x An object of class "\code{bfsl}".
#' @param grid If \code{TRUE} (default) grid lines are plotted.
#' @param ... Further parameters to be passed to the plotting routines.
#'
#' @importFrom graphics abline plot points arrows
#'
#' @export
plot.bfsl = function(x, grid = TRUE, ...)
{
  obj = x
  x = obj$data$x
  y = obj$data$y
  sd_x = obj$data$sd_x
  sd_y = obj$data$sd_y

  plot(x, y, xlim = c(min(x-sd_x), max(x+sd_x)),
       ylim = c(min(y-sd_y), max(y+sd_y)), type = "n", ...)
  if (grid) grid()
  points(x, y, ...)

  # error bars
  arrows(x, y-sd_y, x, y+sd_y, length = 0.05, angle = 90, code = 3)
  arrows(x-sd_x, y, x+sd_x, y, length = 0.05, angle = 90, code = 3)

  # fit line
  abline(coef = obj$coefficients[,1])
}


#' Predict Method for bfsl Model Fits
#'
#' \code{predict.bfsl} predicts future values based on the bfsl fit.
#'
#' @param object Object of class \code{"bfsl"}.
#' @param newdata A data frame with variable \code{x} to predict.
#' If omitted, the fitted values are used.
#' @param interval Type of interval calculation. \code{"none"} or \code{"confidence"}.
#' @param level Confidence level.
#' @param se.fit A switch indicating if standard errors are returned.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return \code{predict.bfsl} produces a vector of predictions or a matrix of
#' predictions and bounds with column names \code{fit}, \code{lwr}, and \code{upr}
#' if interval is set to \code{"confidence"}.
#'
#' If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#' \tabular{ll}{
#' \code{fit} \tab Vector or matrix as above \cr
#' \code{se.fit} \tab Standard error of predicted means
#' }
#'
#' @examples
#' fit = bfsl(pearson_york_data)
#' predict(fit, interval = "confidence")
#' new = data.frame(x = seq(0, 8, 0.5))
#' predict(fit, new, se.fit = TRUE)
#'
#' pred.clim = predict(fit, new, interval = "confidence")
#' matplot(new$x, pred.clim, lty = c(1,2,2), type = "l", xlab = "x", ylab = "y")
#' df = fit$data
#' points(df$x, df$y)
#' arrows(df$x, df$y-df$sd_y, df$x, df$y+df$sd_y,
#'        length = 0.05, angle = 90, code = 3)
#' arrows(df$x-df$sd_x, df$y, df$x+df$sd_x, df$y,
#'        length = 0.05, angle = 90, code = 3)
#'
#' @importFrom stats model.frame model.matrix terms qt
#'
#' @export
predict.bfsl = function(object, newdata, interval = c("none", "confidence"),
                        level = 0.95, se.fit = FALSE, ...) {
  if (missing(newdata) || is.null(newdata)) {
    newdata = data.frame(x = object$data$x)
  }
  if (!is.vector(newdata) || is.list(newdata)) {
    newdata = as.data.frame(newdata)
  }
  if (!("x" %in% colnames(newdata))) {
    stop('No column with name "x" found in newdata.')
  }
  m = model.frame(terms(~x), newdata)
  X = model.matrix(terms(~x), m)
  beta = object$coefficients[,1]
  predictor = drop(X %*% beta)
  interval = match.arg(interval)
  if (se.fit || interval != "none") {
    V = diag(object$coefficients[,2]^2)  # variance covariance matrix
    V[1,2] = object$cov.ab
    V[2,1] = object$cov.ab
    var.fit = rowSums((X %*% V) * X)  # point-wise variance for predicted mean
  }
  if (interval=="confidence") {
    df = object$df.residual  # degrees of freedom
    tfrac = c(-1, 1)*stats::qt((1-level)/2, df, lower.tail = FALSE)  # quantiles of t-distribution
    predictor = cbind(predictor, predictor + outer(sqrt(var.fit), tfrac))
    colnames(predictor) = c("fit", "lwr", "upr")
  }
  if (se.fit) {
    predictor = list(fit = predictor, se.fit = as.numeric(sqrt(var.fit)))
  }
  return(predictor)
}


#' Summary Method for bfsl Results
#'
#' \code{summary} method for class \code{"bfsl"}.
#'
#' @param object An object of class "\code{bfsl}".
#' @param ... Further arguments passed to \code{summary.default}.
#'
#' @export
summary.bfsl = function(object, ...) {

  ans = object
  class(ans) = "summary.bfsl"
  ans
}


#' Print Method for summary.bfsl Objects
#'
#' \code{print} method for class \code{"summary.bfsl"}.
#'
#' @param x An object of class "\code{summary.bfsl}".
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to \code{print.default}.
#'
#' @export
print.summary.bfsl = function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")

  cat("Residuals:\n")
  if(NROW(x$residuals) > 5L) {
    print.default(format(summary(x$residuals), digits = digits), print.gap = 2L,
                  quote = FALSE, ...)
  } else {
    print.default(format(x$residuals, digits = digits), print.gap = 2L,
                  quote = FALSE, ...)
  }

  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE, ...)

  cat("\nGoodness of fit:", format(x$chisq, digits = digits))
  cat("\nChisq-statistic:", format(x$chisq*x$df.residual, digits = digits),
      "on", format(x$df.residual, digits = digits),
      "degrees of freedom")
  cat("\nCovariance of the slope and intercept:", format(x$cov.ab, digits = digits))

  cat("\n\n")
  invisible(x)
}

#' @importFrom generics tidy
#' @export
generics::tidy


#' @importFrom generics glance
#' @export
generics::glance


#' @importFrom generics augment
#' @export
generics::augment


#' Tidy a bfsl Object
#'
#' Broom tidier method to \code{tidy} a bfsl object.
#'
#' @param x A `bfsl` object.
#' @param conf.int Logical indicating whether or not to include
#'   a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level The confidence level to use for the confidence
#'   interval if conf.int = TRUE. Must be strictly greater than 0
#'   and less than 1. Defaults to 0.95, which corresponds to a
#'   95 percent confidence interval.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing component-level
#'   information about the model
#'
#' @examples
#' fit = bfsl(pearson_york_data)
#'
#' tidy(fit)
#'
#' @export
tidy.bfsl = function(x, conf.int = FALSE, conf.level = 0.95, ...) {

  rownames(x$coefficients) = c("(Intercept)", "Slope")
  colnames(x$coefficients) = c("estimate", "std.error")
  result = tibble::as_tibble(x$coefficients, rownames = "term")

  if (conf.int) {
    df = x$df.residual  # degrees of freedom
    tfrac = c(-1, 1)*stats::qt((1-conf.level)/2, df, lower.tail = FALSE)  # quantiles of t-distribution
    ci = x$coefficients[,1] + outer(x$coefficients[,2], tfrac)
    colnames(ci) = c("conf.low", "conf.high")
    ci = tibble::as_tibble(ci, rownames = "term")
    result = dplyr::left_join(result, ci, by = "term")
  }

  result
}


#' Glance at a bfsl Object
#'
#' Broom tidier method to \code{glance} at a bfsl object.
#'
#' @param x A `bfsl` object.
#' @param ... Unused, included for generic consistency only.
#' @return A [tibble::tibble()] with one row and columns:
#' \item{chisq}{The goodness of fit.}
#' \item{p.value}{P-value.}
#' \item{df.residual}{Residual degrees of freedom.}
#' \item{nobs}{Number of observations.}
#' \item{isConv}{Did the fit converge?}
#' \item{iter}{Number of iterations.}
#' \item{finTol}{Final tolerance.}
#'
#' @examples
#' fit = bfsl(pearson_york_data)
#'
#' glance(fit)
#'
#' @export
glance.bfsl = function(x, ...) {
  with(
    summary(x),
    tibble::tibble(
      chisq = chisq,
      df.residual = df.residual,
      nobs = length(data$x),
      isConv = convInfo$isConv,
      iter = convInfo$finIter,
      finTol = convInfo$finTol
    )
  )
}


#' Augment Data with Information from a bfsl Object
#'
#' Broom tidier method to \code{augment} data with information from a bfsl object.
#'
#' @param x A `bfsl` object created by [bfsl::bfsl()]
#' @param data A [base::data.frame()] or [tibble::tibble()] containing all the
#' original predictors used to create x. Defaults to NULL, indicating that
#' nothing has been passed to newdata. If newdata is specified, the data argument
#' will be ignored.
#' @param newdata A [base::data.frame()] or [tibble::tibble()] containing all
#' the original predictors used to create x. Defaults to NULL, indicating that
#' nothing has been passed to newdata. If newdata is specified, the data
#' argument will be ignored.
#' @param ... Unused, included for generic consistency only.
#'
#' @return A [tibble::tibble()] with columns:
#' \item{.fitted}{Fitted or predicted value.}
#' \item{.se.fit}{Standard errors of fitted values.}
#' \item{.resid}{The residuals, that is \code{y} observations minus fitted
#' values. (Only returned if \code{newdata = NULL}).}
#'
#' @examples
#' fit = bfsl(pearson_york_data)
#'
#' augment(fit)
#'
#' @export
augment.bfsl = function(x, data = x$data, newdata = NULL, ...) {
  if (is.null(newdata)) {
    dplyr::bind_cols(tibble::as_tibble(data),
                     tibble::tibble(.fitted = x$fitted.values,
                                    .se.fit = predict.bfsl(x,
                                                      newdata = data,
                                                      se.fit = TRUE)$se.fit,
                                    .resid =  x$residuals))
  } else {
    predictions = predict.bfsl(x, newdata = newdata, se.fit = TRUE)
    dplyr::bind_cols(tibble::as_tibble(newdata),
                     tibble::tibble(.fitted = predictions$fit,
                                    .se.fit = predictions$se.fit))
  }
}

