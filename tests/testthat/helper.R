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
