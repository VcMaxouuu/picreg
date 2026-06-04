# picreg: Variable Selection using the Pivotal Information Criterion

Sparse regression and classification via the Pivotal Information
Criterion (PIC), an alternative to the Bayesian Information Criterion
(BIC), cross-validation, and Lasso-based tuning. The regularization
parameter is selected from a pivotal null-distribution statistic,
eliminating the need for cross-validation and yielding sharper support
recovery. Provides Fast Iterative Shrinkage-Thresholding Algorithm
(FISTA) optimization for the L1, Smoothly Clipped Absolute Deviation
(SCAD), and Minimax Concave Penalty (MCP) penalties across six response
distributions: Gaussian, binomial, Poisson, exponential, Gumbel, and
Cox. Under standard sparsity assumptions, the selector achieves a phase
transition for exact support recovery, analogous to results in
compressed sensing. See Sardy, van Cutsem and van de Geer (2026)
[doi:10.48550/arXiv.2603.04172](https://doi.org/10.48550/arXiv.2603.04172)
.

## See also

Useful links:

- <https://github.com/VcMaxouuu/picreg>

- <https://vcmaxouuu.github.io/picreg/>

- Report bugs at <https://github.com/VcMaxouuu/picreg/issues>

## Author

**Maintainer**: Maxime van Cutsem <maxime.vancutsem@unige.ch>

Authors:

- Maxime van Cutsem <maxime.vancutsem@unige.ch>

- Sylvain Sardy <sylvain.sardy@unige.ch>
