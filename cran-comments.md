## Summary of changes
* Added broom tidier methods `tidy()`, `glance()` and `augment()` for 'bfsl' objects.
* Added `summary.bfsl()` and `print.summary.bfsl()` methods.
* Added `predict.bfsl()` method to calculate confidence intervals (#2).
* Added `bfsl.formula()` method. This allows to plot confidence intervals with `ggplot2`.
* Using `arrows()` instead of `segments()` to plot error bars in `plot.bfsl()`.
* Renamed sample data `pearson_york` to `pearson_york_data`.
* Changed license to MIT.
* Added missing "/n" in `print.bfsl()` (#1).

## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are currently no downstream dependencies for this package.
