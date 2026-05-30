## Resubmission

This is a resubmission of 0.1.1. The previous run on the CRAN Windows
builder reported a WARNING about the rebuilt vignette PDF being larger
than necessary ("compacted 'vignette.pdf' from 455Kb to 155Kb"). In
this version I have switched the vignette output from
`pdf_document` to `rmarkdown::html_vignette`, which eliminates the
PDF compaction issue entirely.

The fixes from the previous resubmission cycle are retained:

* DESCRIPTION: reference written as "Sardy, van Cutsem and van de Geer
  (2026) <doi:10.48550/arXiv.2603.04172>".
* DESCRIPTION: acronyms (PIC, BIC, FISTA, SCAD, MCP) expanded on first
  use.
* cox_partial_log_likelihood.Rd: \value tag added.

## R CMD check results

0 errors | 0 warnings | 1 note

* The spell-check flags acronyms ("FISTA", "MCP", "SCAD", "PIC") and
  proper nouns from the reference. These are intentional.

## Test environments

* local: macOS 15.0 (R 4.6.0)
* win-builder: R-devel and R-release
