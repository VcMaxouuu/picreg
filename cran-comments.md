## Submission

This is an update of the released version 0.1.2.

The CRAN package check reported an "Additional issue (LTO)": a C++ One
Definition Rule violation under link-time optimization. This release
fixes it by giving the file-local helper structs in `src/fista.cpp` and
`src/fista_cox.cpp` internal linkage (anonymous namespaces), so the
`-Wodr` warnings no longer occur. I verified this locally by building
with `-flto -Wodr`.

This release also includes user-facing improvements (see NEWS.md):

* `coef()` now returns a sparse matrix (class "dgCMatrix"), consistent
  with `glmnet`.
* a `summary()` method for fitted models.
* `intercept = FALSE` is now allowed for all families.
* standardized input-validation error messages.

## R CMD check results

0 errors | 0 warnings | 0 notes

The only NOTE seen locally is environmental ("Skipping checking HTML
validation" / "Skipping checking math rendering"), caused by an old
HTML Tidy and the 'V8' package being unavailable on the local machine.
It does not occur on CRAN's check machines.

## Test environments

* local: macOS 15.0 (R 4.6.0), including a build with LTO enabled
  (`-flto -Wodr`) to confirm the reported issue is resolved.
* win-builder: R-devel and R-release
