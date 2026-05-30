## Resubmission

This is a resubmission. In this version I have addressed the points
raised by Benjamin Altmann:

* DESCRIPTION: the reference is now written in the form
  "Sardy, van Cutsem and van de Geer (2026) <doi:10.48550/arXiv.2603.04172>".
* DESCRIPTION: all acronyms (PIC, BIC, FISTA, SCAD, MCP) are now
  expanded on first use.
* cox_partial_log_likelihood.Rd: a \value tag has been added,
  describing the structure and meaning of the returned numeric scalar.

## R CMD check results

0 errors | 0 warnings | 1 note

* The spell-check flags acronyms ("FISTA", "MCP", "SCAD", "PIC") and
  proper nouns from the reference. These are intentional.

## Test environments

* local: macOS 15.0 (R 4.6.0)
* win-builder: R-devel and R-release
