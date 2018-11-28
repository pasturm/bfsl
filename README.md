# R package bfst: Best-fit Straight Line
[![Travis build status](https://travis-ci.org/pasturm/bfsl.svg?branch=master)](https://travis-ci.org/pasturm/bfsl)

### How to fit a straight line through a set of points with errors in both coordinates?

The solution for the best-fit straight line to independent points with normally distributed errors in both _x_ and _y_ is known from York (1966, 1969, 2004). It provides unbiased estimates of the intercept, slope and standard errors of the best-fit straight line, even when the _x_ and _y_ errors are correlated.

Surprisingly, as Wehr and Saleska (2017) point out, York's solution has escaped the attention of many statisticians and other scientists that are writing on straight-line fitting with errors in both _x_ and _y_ (also known as Model II regressions, errors-in-variables models or measurement error models).

Other commonly used least-squares estimation methods, such as ordinary least-squares regression, orthogonal distance regression (also known as major axis regression), and geometric mean regression (also known as reduced major axis or standardised major axis regression), are all special cases of York’s solution and only valid under particular measurement conditions.

The bfsl package implements York's general solution and provides the best-fit straight line of bivariate data with errors in both coordinates.

### Installation
```
if (!require("devtools")) { install.packages("devtools") }
devtools::install_github("pasturm/bfsl")
```

### Release notes
See the [NEWS file](https://github.com/pasturm/bfsl/blob/master/NEWS.md) for latest release notes.

### References

York, D. (1966). Least-squares fitting of a straight line. _Canadian Journal of Physics_, 44(5), 1079–1086, https://doi.org/10.1139/p66-090

York, D. (1969). Least squares fitting of a straight line with correlated errors. _Earth and Planetary Science Letters_, 5, 320–324,
https://doi.org/10.1016/S0012-821X(68)80059-7

York D. and N.M. Evensen (2004). Unified equations for the slope, intercept, and 
standard errors of the best straight line, _American Journal of Physics_, 72, 367-375, https://doi.org/10.1119/1.1632486

Wehr, R. and Saleska, S. R. (2017) The long-solved problem of the best-fit straight line: application to isotopic mixing lines, _Biogeosciences_, 14, 17-29, https://doi.org/10.5194/bg-14-17-2017
