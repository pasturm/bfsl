context("test-bfsl")

test_that("bfsl works with pearso_york sample data", {
  x = pearson_york$x
  y = pearson_york$y
  sd_x = 1/sqrt(pearson_york$w_x)
  sd_y = 1/sqrt(pearson_york$w_y)
  fit = bfsl(x, y, sd_x, sd_y)
  expect_equal(fit$convInfo$isConv, TRUE)
  expect_equal(fit$coefficients[1,1], 5.479910224)
  expect_equal(fit$coefficients[2,1], -0.480533407)
  expect_equal(fit$coefficients[1,2], 0.294970735)
  expect_equal(fit$coefficients[2,2], 0.0579850090)
  expect_equal(fit$chisq, 1.483294149)
})

test_that("bfsl gives ordinary least squares solution", {
  fit = bfsl(pearson_york$x, pearson_york$y)
  ols = summary(lm(y ~ x, pearson_york))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
  # weighted ols
  fit = bfsl(pearson_york$x, pearson_york$y, 0, 1/sqrt(pearson_york$w_y))
  ols = summary(lm(y ~ x, pearson_york, weights = w_y))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
})

test_that("bfsl gives geometric mean regression solution", {
  fit = bfsl(pearson_york$x, pearson_york$y, sd(pearson_york$x), sd(pearson_york$y))
  ols = summary(lm(y ~ x, pearson_york))
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1]/sqrt(ols$r.squared))
  expect_equal(fit$coefficients[2,1], -sd(pearson_york$y)/sd(pearson_york$x))
})

test_that("bfsl gives Deming regression solution", {
  sd_x = 3
  sd_y = 7
  fit = bfsl(pearson_york$x, pearson_york$y, sd_x, sd_y)
  Sxx = sum((pearson_york$x-mean(pearson_york$x))^2)
  Syy = sum((pearson_york$y-mean(pearson_york$y))^2)
  Sxy = sum((pearson_york$x-mean(pearson_york$x))*(pearson_york$y-mean(pearson_york$y)))
  lambda = sd_y^2/sd_x^2
  b_deming = (Syy - lambda*Sxx + sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy)
  expect_equal(fit$coefficients[2,1], b_deming)
})

test_that("bfsl works with data frames, lists, matrices", {
  x = pearson_york$x
  y = pearson_york$y
  sd_x = 1/sqrt(pearson_york$w_x)
  sd_y = 1/sqrt(pearson_york$w_y)
  fit1 = bfsl(x, y, sd_x, sd_y)
  fit2 = bfsl(pearson_york)
  fit3 = bfsl(as.matrix(pearson_york))
  fit4 = bfsl(as.array(as.matrix(pearson_york)))
  fit5 = bfsl(list(x = x, y = y, sd_x = sd_x, sd_y = sd_y))
  fit6 = bfsl(as.data.frame(x), y, sd_x, sd_y)
  expect_equal(fit2$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit3$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit4$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit5$coefficients[1,1], fit1$coefficients[1,1])
  expect_equal(fit6$coefficients[1,1], fit1$coefficients[1,1])
})
