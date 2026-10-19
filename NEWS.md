# scramble 0.2.0

## Changes to numerical results

Version 0.1.0 reproduces the results of the code used in published papers
(up to ~1e-10). The following fixes change results; each is a separate commit
so that results can be traced to the version that produced them.

* `loss_type = "LTS"`: the trimming threshold was the h-quantile of the
  *signed* residuals, compared with `<` against the *squared* residuals. The
  number of kept cells therefore depended on the scale of the residuals
  (about 58% instead of 75% for standard normal residuals with `param = 0.75`,
  fewer for larger residuals, and none for `param < 0.5`). The threshold is now
  the h-quantile of the squared residuals of each variable and cells with
  squared residual `<=` threshold are kept, so that the proportion `param` of
  cells per variable is used as documented. Results with `loss_type = "LTS"`
  change; other loss types are unaffected.

* `transformation = "wrapping"`: the data, already robustly standardized
  with `center = TRUE` / `scale = TRUE`, were passed to `cellWise::wrap()`
  together with the location and scale of the *raw* data, so they were
  standardized a second time. Depending on the location and scale of the raw
  variables this distorted the initial estimate or produced missing values
  (and an error in `svd()`). `wrap()` now gets location 0 and scale 1 for
  standardized data, and otherwise the median and Qn of the data it is
  applied to. In simulations (n = 100, p = 20, 5% cellwise outliers) the
  median angle between the initial and the true loadings dropped from
  59/81 to 13/20 degrees. Results with `transformation = "wrapping"` change,
  except with `center = FALSE, scale = FALSE` without groups.

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
