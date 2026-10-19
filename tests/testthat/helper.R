skip_if_no_torch <- function() {
  skip_if_not(torch::torch_is_installed(), "torch backend (libtorch) not installed")
}

# small data set with a two-dimensional structure and a few cellwise outliers
make_data <- function(n = 40, p = 6, n_outliers = 10) {
  set.seed(1)
  x <- matrix(rnorm(n * p), n, p)
  x[sample(length(x), n_outliers)] <- 15
  x
}

# pcaSCRAMBLE with fixed penalty and few epochs; convergence warnings expected
fit_quick <- function(x = make_data(), ..., epochs = 150) {
  suppressWarnings(
    pcaSCRAMBLE(x, k = 2, alpha_x = 0.5, epochs = epochs,
                penalties = list(pen_x = 0.1), ...)
  )
}

is_orthonormal <- function(v, tol = 1e-6) {
  isTRUE(all.equal(crossprod(v), diag(ncol(v)), tolerance = tol,
                   check.attributes = FALSE))
}

# evaluate expr, muffling and returning all warning messages
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

# configurations of the reference-value test; the reference values are created
# by fixtures/make-reference.R and must only be regenerated deliberately
reference_configs <- list(
  L2 = list(),
  Huber = list(loss_type = "Huber"),
  wrapping_LTS = list(transformation = "wrapping", loss_type = "LTS",
                      param = 0.75),
  spearman_noscale = list(transformation = "spearman", scale = FALSE),
  groups_L2 = list(groups = rep(1:2, 20)),
  groups_LTS = list(groups = rep(1:2, 20), loss_type = "LTS")
)
reference_fields <- c("a", "a_ortho", "measure", "measure_ortho", "phi")
