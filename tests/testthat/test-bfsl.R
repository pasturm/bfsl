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

test_that("bfsl gives ols solution", {
  fit = bfsl(pearson$x, pearson$y)
  ols = summary(lm(pearson$y~pearson$x))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
  # weighted ols
  fit = bfsl(pearson$x, pearson$y, 0, pearson$sd_y)
  ols = summary(lm(pearson$y~pearson$x, weights = 1/pearson$sd_y^2))
  expect_equal(fit$coefficients[1,1], ols$coefficients[1,1])
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1])
  expect_equal(fit$coefficients[1,2]*sqrt(fit$chisq), ols$coefficients[1,2])
  expect_equal(fit$coefficients[2,2]*sqrt(fit$chisq), ols$coefficients[2,2])
})

test_that("bfsl gives gmr solution", {
  fit = bfsl(pearson$x, pearson$y, sd(pearson$x), sd(pearson$y))
  ols = summary(lm(pearson$y~pearson$x))
  expect_equal(fit$coefficients[2,1], ols$coefficients[2,1]/sqrt(ols$r.squared))
  expect_equal(fit$coefficients[2,1], -sd(pearson$y)/sd(pearson$x))
})

test_that("bfsl gives odr solution", {
  fit = bfsl(pearson$x, pearson$y, 1, 1)
  Sxx = sum((pearson$x-mean(pearson$x))^2)
  Syy = sum((pearson$y-mean(pearson$y))^2)
  Sxy = sum((pearson$x-mean(pearson$x))*(pearson$y-mean(pearson$y)))
  # b_ols = Sxy/Sxx
  b_odr = (Syy - Sxx + sqrt((Syy-Sxx)^2+4*Sxy^2))/(2*Sxy)
  expect_equal(fit$coefficients[2,1], b_odr)
})
