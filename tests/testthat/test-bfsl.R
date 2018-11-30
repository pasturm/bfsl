context("test-bfsl")

test_that("bfsl works with pearson sample data", {
  fit = bfsl(pearson$x, pearson$y, pearson$sd_x, pearson$sd_y)
  expect_equal(fit$convInfo$isConv, TRUE)
  expect_equal(fit$coefficients[1,1], 5.475612048)
  expect_equal(fit$coefficients[2,1], -0.479411379)
  expect_equal(fit$coefficients[1,2], 0.2927733369)
  expect_equal(fit$coefficients[2,2], 0.0576051364)
  expect_equal(fit$chisq, 1.507862028)
})

test_that("bfsl gives ordinary least squares solution", {
  fit = bfsl(pearson$x, pearson$y)
  ols = summary(lm(y ~ x, pearson))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
  # weighted ols
  fit = bfsl(pearson$x, pearson$y, 0, pearson$sd_y)
  ols = summary(lm(y ~ x, pearson, weights = 1/sd_y^2))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
})

test_that("bfsl gives geometric mean regression solution", {
  fit = bfsl(pearson$x, pearson$y, sd(pearson$x), sd(pearson$y))
  ols = summary(lm(y ~ x, pearson))
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1]/sqrt(ols$r.squared))
  expect_equal(fit$coefficients[2,1], -sd(pearson$y)/sd(pearson$x))
})

test_that("bfsl gives Deming regression solution", {
  sd_x = 3
  sd_y = 7
  fit = bfsl(pearson$x, pearson$y, sd_x, sd_y)
  Sxx = sum((pearson$x-mean(pearson$x))^2)
  Syy = sum((pearson$y-mean(pearson$y))^2)
  Sxy = sum((pearson$x-mean(pearson$x))*(pearson$y-mean(pearson$y)))
  lambda = sd_y^2/sd_x^2
  b_deming = (Syy - lambda*Sxx + sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy)
  expect_equal(fit$coefficients[2,1], b_deming)
})

