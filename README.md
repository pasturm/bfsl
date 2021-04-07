
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package bfsl: Best-fit Straight Line

<!-- badges: start -->

[![R build
status](https://github.com/pasturm/bfsl/workflows/R-CMD-check/badge.svg)](https://github.com/pasturm/bfsl/actions)
[![codecov](https://codecov.io/gh/pasturm/bfsl/branch/master/graph/badge.svg)](https://codecov.io/gh/pasturm/bfsl)
[![CRAN
version](https://www.r-pkg.org/badges/version-last-release/bfsl)](https://cran.r-project.org/package=bfsl)
[![CRAN
downloads/month](https://cranlogs.r-pkg.org/badges/bfsl)](https://cran.r-project.org/package=bfsl)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/bfsl)](https://cran.r-project.org/package=bfsl)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

### How to fit a straight line through a set of points with errors in both coordinates?

The solution for the best-fit straight line to independent points with
normally distributed errors in both *x* and *y* is known e.g. from York
([1966](#york66), [1968](#york68), [2004](#york04)). It provides
unbiased estimates of the intercept, slope and standard errors of the
best-fit straight line, even when the *x* and *y* errors are correlated.

The bfsl package implements York’s general solution and provides the
best-fit straight line of bivariate data with errors in both
coordinates.

Other commonly used least-squares estimation methods, such as ordinary
least-squares regression, orthogonal distance regression (also called
major axis regression), geometric mean regression (also called reduced
major axis or standardised major axis regression) or Deming regression
are all special cases of York’s solution and only valid under particular
measurement conditions.

## Example

<!-- # ```{r print} -->
<!-- # library(bfsl) -->
<!-- # bfsl(pearson_york_data) -->
<!-- # ``` -->

``` r
library(bfsl)
fit = bfsl(pearson_york_data)
print(fit)
#> Best-fit straight line
#> 
#>            Estimate  Std. Error
#> Intercept   5.47991   0.29497  
#> Slope      -0.48053   0.05799  
#> 
#> Goodness of fit:
#> 1.483

plot(fit)
ols = bfsl(pearson_york_data, sd_x = 0, sd_y = 1)
abline(coef = ols$coef[,1], lty = 2)
legend("topright", c("ordinary least squares", "best-fit straight line"), lty = c(2,1))
```

<img src="man/figures/README-plot-1.png" width="75%" />

Confindence interval:

``` r
df = as.data.frame(fit$data)
newx = seq(min(df$x-df$sd_x), max(df$x+df$sd_x), length.out = 100)
preds = predict(fit, newdata = data.frame(x=newx), interval = 'confidence')

# plot
plot(y ~ x, data = df, type = 'n', xlim = c(min(df$x-df$sd_x), max(df$x+df$sd_x)),
     ylim = c(min(df$y-df$sd_y), max(df$y+df$sd_y)), las = 1)
grid()
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey90', border = NA)
abline(coef = fit$coef[,1], lty = 1)
points(df$x, df$y)
arrows(df$x, df$y-df$sd_y, df$x, df$y+df$sd_y, length = 0.05, angle = 90, code = 3)
arrows(df$x-df$sd_x, df$y, df$x+df$sd_x, df$y, length = 0.05, angle = 90, code = 3)
```

<img src="man/figures/README-confidence-1.png" width="75%" />

``` r
# with ggplot2
library(ggplot2)
ggplot(data = df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = bfsl) +
  geom_errorbar(aes(ymin = y-sd_y, ymax = y+sd_y), width = 0.05) +
  geom_errorbarh(aes(xmin = x-sd_x, xmax = x+sd_x), height = 0.05)
#> `geom_smooth()` using formula 'y ~ x'
```

<img src="man/figures/README-ggplot-1.png" width="60%" />

## Installation

Install bfsl from CRAN with:

``` r
install.packages("bfsl")
```

Install the development version from GitHub with:

``` r
if (!require("remotes")) { install.packages("remotes") }
remotes::install_github("pasturm/bfsl")
```

See the [NEWS file](https://github.com/pasturm/bfsl/blob/master/NEWS.md)
for latest release notes.

## References

<a name="york66"></a>York, D. (1966). Least-squares fitting of a
straight line. *Canadian Journal of Physics*, 44(5), 1079–1086,
<https://doi.org/10.1139/p66-090>

<a name="york68"></a>York, D. (1968). Least squares fitting of a
straight line with correlated errors. *Earth and Planetary Science
Letters*, 5, 320–324, <https://doi.org/10.1016/S0012-821X(68)80059-7>

<a name="williamson68"></a>Williamson, J. H. (1968). Least-squares
fitting of a straight line, *Canadian Journal of Physics*, 46(16),
1845-1847, <https://doi.org/10.1139/p68-523>

<a name="york04"></a>York, D. et al. (2004). Unified equations for the
slope, intercept, and standard errors of the best straight line,
*American Journal of Physics*, 72, 367-375,
<https://doi.org/10.1119/1.1632486>

<a name="cantrell08"></a>Cantrell, C. A. (2008). Technical Note: Review
of methods for linear least-squares fitting of data and application to
atmospheric chemistry problems, *Atmospheric Chemistry and Physics*, 8,
5477-5487, <https://www.atmos-chem-phys.net/8/5477/2008/>

<a name="wehr17"></a>Wehr, R. and Saleska, S. R. (2017). The long-solved
problem of the best-fit straight line: application to isotopic mixing
lines, *Biogeosciences*, 14, 17-29,
<https://doi.org/10.5194/bg-14-17-2017>
