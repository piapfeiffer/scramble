# scramble

**Cellwise robust and sparse principal component analysis in R**

`scramble` implements the sparse and cellwise robust PCA method proposed in

> Pfeiffer, P., Vana-Gür, L., & Filzmoser, P. (2025). Cellwise robust and
> sparse principal component analysis. *Advances in Data Analysis and
> Classification*. <https://doi.org/10.1007/s11634-025-00656-3>

Classical PCA, and most robust variants, treat whole observations as either
clean or outlying. In high-dimensional data, contamination often affects
single cells of the data matrix instead, and these cellwise outliers spread
across most rows. `scramble` handles this by

- replacing the squared reconstruction error with a **robust loss** (Huber,
  Tukey or least trimmed squares),
- adding an **elastic net penalty** for sparse, interpretable loadings,
- starting from a **robust initial estimate** based on a rank-based
  (Spearman) or wrapping transformation, and
- optimizing directly on the **Stiefel manifold** with Riemannian gradient
  descent, implemented with [`torch`](https://torch.mlverse.org/) for R.

The penalty parameter can be fixed or tuned automatically by Bayesian
optimization.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("piapfeiffer/scramble")
```

`scramble` uses `torch`, which needs to download its backend (libtorch) once
after installation:

```r
torch::install_torch()
```

## Example

```r
library(scramble)

set.seed(1)
x <- matrix(rnorm(100 * 10), 100, 10)
x[, 2] <- x[, 1] + rnorm(100, sd = 0.2)
x[, 3] <- x[, 1] + rnorm(100, sd = 0.2)
x[sample(length(x), 30)] <- 10   # cellwise outliers

res <- pcaSCRAMBLE(x,
                   transformation = "wrapping",  # robust starting value
                   loss_type = "Huber",          # robust loss
                   alpha_x = 0.5,                # elastic net mixing
                   k = 2,                        # number of components
                   penalties = list(pen_x = 0.5))

round(res$a, 2)   # sparse loadings: the first component picks out variables 1-3
plot(res)
```

Leave out `penalties` to tune the penalty by Bayesian optimization. See
`?pcaSCRAMBLE` and `vignette("demo", package = "scramble")` for all options.

## Citation

If you find this package useful in your work, please cite the paper:

> Pfeiffer, P., Vana-Gür, L., & Filzmoser, P. (2025). Cellwise robust and
> sparse principal component analysis. *Advances in Data Analysis and
> Classification*. <https://doi.org/10.1007/s11634-025-00656-3>

```bibtex
@article{pfeiffer2025cellwise,
  title   = {Cellwise robust and sparse principal component analysis},
  author  = {Pfeiffer, Pia and Vana-G{\"u}r, Laura and Filzmoser, Peter},
  journal = {Advances in Data Analysis and Classification},
  year    = {2025},
  doi     = {10.1007/s11634-025-00656-3}
}
```

You can also get the reference from R with `citation("scramble")`.

## License

MIT © Pia Pfeiffer. The Stiefel manifold optimizer is adapted from the
optimizer template of [`torch`](https://torch.mlverse.org/) (MIT, © Daniel
Falbel and RStudio); see `inst/COPYRIGHTS`.
