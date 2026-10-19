# scramble 0.1.0

* First CRAN release. The user-facing API (`pcaSCRAMBLE()`, `train_PCA_SVD()`,
  `principal_angles()`, `TPR()`, `TNR()`, `plot()`) is unchanged.
* `pcaSCRAMBLE()` validates its inputs and gives informative errors for
  missing or invalid `alpha_x`, unknown `transformation` / `loss_type`, and
  `k > min(n, p)` (previously cryptic errors).
* `loss_type = "LTS"` with the default `param = NA` now uses `param = 0.75`
  instead of failing.
* If Bayesian optimization of the penalty fails, the upper bound of `bounds`
  is used as documented (previously `NA`, which made training fail), and a
  warning is raised instead of printing.
* `plot()` on a `PCA_result` no longer fails, and restores the graphical
  parameters afterwards.
* The `loss` trace now records the first epoch.
* Removed the dependency on 'rstackdeque'.
* Training is faster: the data tensor is created once instead of in every
  epoch, and the robust variances are computed once after training instead of
  at every improvement. Results are bit-identical to before.
* Training is up to ~5x faster for many variables: the Cayley step of the
  Stiefel optimizer solves a (2k x 2k) system via the Woodbury identity
  instead of a (p x p) system, the reconstruction is computed as
  `(X V) V^T`, and the elastic net penalty is computed for all components
  at once. These are mathematically identical, but floating-point rounding
  differs: loadings change by about 1e-10 (sparsity patterns and selected
  penalties were unchanged in our comparisons).
