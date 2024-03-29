# bfsl 0.2.1

* Removed p-value of chi-square statistic in `summary.bfsl()` and `glance.bfsl()`.
* Added t-statistic and corresponding p-value of the fitted coefficients in
  `summary.bfsl()`.

# bfsl 0.2.0

* Added `summary.bfsl()` and `print.summary.bfsl()` methods.
* Added broom tidier methods `tidy()`, `glance()` and `augment()` for bfsl objects.
* Updated `print.bfsl()`.

# bfsl 0.1.1

* Updated README.md.
* Using `arrows()` instead of `segments()` to plot error bars in `plot.bfsl()`.
* Renamed sample data `pearson_york` to `pearson_york_data`.
* Changed license to MIT.
* Added missing "/n" in `print.bfsl()` (#1).
* Added `predict.bfsl()` method to calculate confidence intervals (#2).
* Added `bfsl.formula()` method. This allows to plot confidence intervals
  with `ggplot2`.

# bfsl 0.1.0

* Initial release.
