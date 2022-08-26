context("test-bfsl")

test_that("bfsl works with pearson_york sample data", {
  x = pearson_york_data$x
  y = pearson_york_data$y
  sd_x = 1/sqrt(pearson_york_data$w_x)
  sd_y = 1/sqrt(pearson_york_data$w_y)
  fit = bfsl(x, y, sd_x, sd_y)
  expect_equal(fit$convInfo$isConv, TRUE)
  expect_equal(fit$coefficients[1,1], 5.479910224)
  expect_equal(fit$coefficients[2,1], -0.480533407)
  expect_equal(fit$coefficients[1,2], 0.294970735)
  expect_equal(fit$coefficients[2,2], 0.0579850090)
  expect_equal(fit$chisq, 1.483294149)
})

test_that("bfsl gives ordinary least squares solution", {
  fit = bfsl(pearson_york_data$x, pearson_york_data$y)
  ols = summary(lm(y ~ x, pearson_york_data))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
  # weighted ols
  fit = bfsl(pearson_york_data$x, pearson_york_data$y, 0, 1/sqrt(pearson_york_data$w_y))
  ols = summary(lm(y ~ x, pearson_york_data, weights = w_y))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
})

test_that("bfsl gives geometric mean regression solution", {
  fit = bfsl(pearson_york_data$x, pearson_york_data$y, sd(pearson_york_data$x), sd(pearson_york_data$y))
  ols = summary(lm(y ~ x, pearson_york_data))
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1]/sqrt(ols$r.squared))
  expect_equal(fit$coefficients[2,1], -sd(pearson_york_data$y)/sd(pearson_york_data$x))
})

test_that("bfsl gives Deming regression solution", {
  sd_x = 3
  sd_y = 7
  fit = bfsl(pearson_york_data$x, pearson_york_data$y, sd_x, sd_y)
  Sxx = sum((pearson_york_data$x-mean(pearson_york_data$x))^2)
  Syy = sum((pearson_york_data$y-mean(pearson_york_data$y))^2)
  Sxy = sum((pearson_york_data$x-mean(pearson_york_data$x))*(pearson_york_data$y-mean(pearson_york_data$y)))
  lambda = sd_y^2/sd_x^2
  b_deming = (Syy - lambda*Sxx + sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy)
  expect_equal(fit$coefficients[2,1], b_deming)
})

test_that("bfsl works with data frames, lists, matrices, formulas", {
  x = pearson_york_data$x
  y = pearson_york_data$y
  sd_x = 1/sqrt(pearson_york_data$w_x)
  sd_y = 1/sqrt(pearson_york_data$w_y)
  fit1 = bfsl(x, y, sd_x, sd_y)
  fit2 = bfsl(pearson_york_data)
  fit3 = bfsl(as.matrix(pearson_york_data))
  fit4 = bfsl(as.array(as.matrix(pearson_york_data)))
  fit5 = bfsl(list(x = x, y = y, sd_x = sd_x, sd_y = sd_y))
  fit6 = bfsl(as.data.frame(x), y, sd_x, sd_y)
  fit7 = bfsl(y ~ x, data = pearson_york_data)
  expect_equal(fit2$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit3$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit4$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit5$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit6$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit7$coefficients[1,1], fit1$coefficients[1,1])
})

test_that("print method works", {
  expect_output(print(bfsl(pearson_york_data)),
"
Call:
bfsl.default(x = pearson_york_data)

Coefficients:\n             Estimate  Std. Error
(Intercept)   5.47991   0.29497  \nSlope        -0.48053   0.05799  \n"
, fixed=TRUE)
})

test_that("plot method does not create an error", {
  expect_equal(plot(bfsl(pearson_york_data)), NULL)
})

test_that("predict method works", {
  fit = bfsl(pearson_york_data)
  pred = predict(fit, interval = 'confidence', se.fit = TRUE)
  expect_equal(as.numeric(pred$fit[1,]), c(5.47991022, 4.79970649, 6.16011396))
  expect_equal(pred$se.fit[1], 0.29497074)
  fit = bfsl(pearson_york_data$x, pearson_york_data$y, sd_x = 0, sd_y = 1)
  pred = predict(fit, interval = 'confidence', se.fit = TRUE)
  ols = lm(y ~ x, pearson_york_data)
  pred_ols = predict(ols, interval = 'confidence', se.fit = TRUE)
  expect_equal(pred$se.fit*sqrt(fit$chisq), pred_ols$se.fit)
})

test_that("summary method works", {
  fit = bfsl(pearson_york_data)
  s = summary(fit)
  expect_output(print(summary(fit)),
"
Call:
bfsl.default(x = pearson_york_data)

Residuals:\n    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  \n-0.42396  -0.20650   0.14745   0.05573   0.34804   0.42009  \n
Coefficients:\n            Estimate Std. Error t value Pr(>|t|)    \n(Intercept)  5.47991    0.29497  18.578 7.27e-08 ***\nSlope       -0.48053    0.05799  -8.287 3.38e-05 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Goodness of fit: 1.483
Chisq-statistic: 11.87 on 8 degrees of freedom
Covariance of the slope and intercept: -0.01651
"
, fixed=TRUE)
})

test_that("tidy method works" , {
  fit = bfsl(pearson_york_data)
  tmp = tidy(fit)
  expect_equal(tmp$term[1], "(Intercept)")
  tmp = tidy(fit, conf.int = TRUE)
  expect_equal(tmp$conf.high[2], -0.346819737)
  tmp = tidy(fit, conf.int = TRUE, conf.level = 0.99)
  expect_equal(tmp$conf.high[2], -0.285971243)
  tmp = bfsl(pearson_york_data$x, pearson_york_data$y, sd_x = 0, sd_y = 1)
  fit = bfsl(pearson_york_data$x, pearson_york_data$y, sd_x = 0, sd_y = sqrt(tmp$chisq))
  tmp = tidy(fit, conf.int = TRUE)
  ols = lm(y ~ x, pearson_york_data)
  conf.low.ols = confint(ols)[,1]
  names(conf.low.ols) = NULL
  conf.low.bfsl = tmp$conf.low
  expect_equal(conf.low.bfsl, conf.low.ols)
})

test_that("glance method works" , {
  fit = bfsl(pearson_york_data)
  tmp = glance(fit)
  expect_equal(tmp$chisq, 1.483294149232)
})

test_that("augment method works" , {
  fit = bfsl(pearson_york_data)
  tmp = augment(fit)
  expect_equivalent(tmp$.resid, fit$residuals)
  expect_equivalent(tmp$.fitted, fit$fitted.values)
  newdata = data.frame(x = c(2:7))
  tmp = augment(fit, newdata = newdata)
  pred = predict(fit, newdata, se.fit = TRUE)
  expect_equivalent(tmp$.fitted, pred$fit)
  expect_equivalent(tmp$.se.fit, pred$se.fit)
})
