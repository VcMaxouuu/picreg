## Submission

This is an update of the released version 0.1.3.

### Changes in this version

* `predict()` for Cox fits now accepts `type = "survival"`, returning
  subject-specific survival curves. This replaces the standalone
  `predict_survival_function()`, which has been removed.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS 15.0 (R 4.6.0)
* win-builder: R-devel and R-release
