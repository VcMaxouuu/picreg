# Articles

A selection of published work on the Pivotal Information Criterion (PIC)
and the methods implemented in `picreg`.

## The Pivotal Information Criterion

Sardy, van Cutsem & van de Geer (2026).
[arXiv:2603.04172](https://arxiv.org/abs/2603.04172)

The companion paper that introduces PIC. It develops the pivotal
detection-boundary principle and the (asymptotic) pivotality of the null
gradient statistic, then studies support recovery and the
phase-transition behavior empirically, benchmarking PIC against
competing methods on real-world datasets.

## A pivotal transform for the high-dimensional location-scale model

van de Geer, Sardy & van Cutsem (2025).
[arXiv:2512.18705](https://arxiv.org/abs/2512.18705)

The theoretical backbone of PIC in the location-scale setting: the
exponential transformation of the log-likelihood whose \ell_1-penalty
tuning parameter does not depend on the unknown scale parameter,
generalizing the square-root Lasso beyond quadratic loss. The tuning
parameter can asymptotically be taken at the detection edge, and the
paper proves an oracle inequality, variable-selection consistency, and
asymptotic efficiency of the scale and intercept estimators, with
examples including the Subbotin and Gumbel distributions.
