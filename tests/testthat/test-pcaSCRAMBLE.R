test_that("pcaSCRAMBLE returns a PCA_result with the documented structure", {
  skip_if_no_torch()
  x <- make_data()
  res <- fit_quick(x)
  expect_s3_class(res, "PCA_result")
  expect_named(res, c("a", "a_ortho", "u", "s", "loss", "measure",
                      "measure_ortho", "explained_var", "phi", "pen_x",
                      "summary", "wt", "robscale", "robcenter", "alpha_x"))
  expect_equal(dim(res$a), c(6, 2))
  expect_equal(dim(res$a_ortho), c(6, 2))
  expect_equal(dim(res$phi), c(40, 2))
  expect_length(res$measure, 2)
  expect_false(is.unsorted(rev(res$measure)))
  expect_equal(res$explained_var[2], 1)
  expect_equal(res$pen_x, 0.1)
  expect_equal(res$alpha_x, 0.5)
  expect_true(is_orthonormal(res$a_ortho))
  expect_equal(colSums(res$a^2), c(1, 1), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("scores are the standardized data times the loadings", {
  skip_if_no_torch()
  x <- make_data()
  res <- fit_quick(x)
  z <- scale(x, center = res$robcenter, scale = res$robscale)
  expect_equal(res$phi, z %*% res$a, ignore_attr = TRUE)
})

test_that("the loss trace starts at the first epoch", {
  skip_if_no_torch()
  res <- fit_quick(epochs = 20)
  expect_length(res$loss, 20)
  expect_false(anyNA(res$loss))
})

test_that("all transformations and loss types run", {
  skip_if_no_torch()
  x <- make_data()
  for (transformation in c("identity", "spearman", "wrapping")) {
    for (loss_type in c("L2", "Huber", "Tukey", "LTS")) {
      res <- fit_quick(x, transformation = transformation,
                       loss_type = loss_type, epochs = 30)
      expect_s3_class(res, "PCA_result")
      expect_false(anyNA(res$a), label = paste(transformation, loss_type))
    }
  }
})

test_that("LTS uses param = 0.75 by default", {
  skip_if_no_torch()
  x <- make_data()
  expect_equal(fit_quick(x, loss_type = "LTS")$a,
               fit_quick(x, loss_type = "LTS", param = 0.75)$a)
})

test_that("groups, uncentered and unscaled data are supported", {
  skip_if_no_torch()
  x <- make_data()
  expect_s3_class(fit_quick(x, groups = rep(1:2, 20), epochs = 30), "PCA_result")
  expect_s3_class(fit_quick(x, groups = rep(1:2, 20), loss_type = "LTS",
                            epochs = 30), "PCA_result")
  expect_s3_class(fit_quick(x, center = FALSE, scale = FALSE, epochs = 30),
                  "PCA_result")
})

test_that("data frames give the same result as matrices", {
  skip_if_no_torch()
  x <- make_data()
  expect_equal(fit_quick(as.data.frame(x), epochs = 30)$a,
               fit_quick(x, epochs = 30)$a, ignore_attr = TRUE)
})

test_that("k defaults to min(n, p) and is capped with a warning", {
  skip_if_no_torch()
  x <- make_data(n = 5, p = 8, n_outliers = 0)
  res <- suppressWarnings(pcaSCRAMBLE(x, alpha_x = 0.5, epochs = 20,
                                      penalties = list(pen_x = 0.1)))
  expect_equal(ncol(res$a), 5)
  out <- collect_warnings(pcaSCRAMBLE(x, k = 7, alpha_x = 0.5, epochs = 20,
                                      penalties = list(pen_x = 0.1)))
  expect_true(any(grepl("k has to be <= min", out$warnings)))
  expect_equal(ncol(out$value$a), 5)
})

test_that("k = 1 works", {
  skip_if_no_torch()
  res <- suppressWarnings(pcaSCRAMBLE(make_data(), k = 1, alpha_x = 0.5,
                                      epochs = 30, penalties = list(pen_x = 0.1)))
  expect_length(res$measure, 1)
  expect_equal(dim(res$phi), c(40, 1))
})

test_that("invalid input gives informative errors", {
  x <- make_data()
  expect_error(pcaSCRAMBLE(x), "alpha_x")
  expect_error(pcaSCRAMBLE(x, alpha_x = 2), "alpha_x")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, transformation = "foo"),
               "transformation")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, loss_type = "foo"), "loss_type")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, loss_type = "LTS", param = 2),
               "param")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, groups = 1:3), "groups")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, k = 0), "k")
  expect_error(pcaSCRAMBLE(x, alpha_x = 0.5, bounds = 1), "bounds")
  expect_error(pcaSCRAMBLE(matrix(letters[1:4], 2), alpha_x = 0.5), "numeric")
})

test_that("plot works and restores graphical parameters", {
  skip_if_no_torch()
  res <- fit_quick(epochs = 30)
  pdf(NULL)
  on.exit(dev.off())
  mfrow <- par("mfrow")
  expect_invisible(plot(res))
  expect_equal(par("mfrow"), mfrow)
})

test_that("results match the reference values", {
  # guards against unintended numerical changes; see
  # fixtures/make-reference.R for how (and when) to regenerate the values
  skip_if_no_torch()
  skip_on_cran()
  x <- make_data()
  ref <- readRDS(test_path("fixtures", "reference.rds"))
  expect_setequal(names(ref), names(reference_configs))
  for (cfg in names(ref)) {
    res <- do.call(fit_quick, c(list(x), reference_configs[[cfg]]))
    expect_equal(res[reference_fields], ref[[cfg]], tolerance = 1e-6,
                 label = cfg)
  }
})

test_that("wrapping is invariant to the location and scale of the variables", {
  skip_if_no_torch()
  x <- make_data()
  y <- sweep(x, 2, c(1, 2, 5, 0.1, 3, 10), `*`) + 7
  fit <- function(x, ...) fit_quick(x, transformation = "wrapping", epochs = 30, ...)
  expect_equal(fit(y)$a, fit(x)$a, tolerance = 1e-6)
  # only scaled: invariant to scaling (location matters, the data are not centered)
  expect_equal(fit(x * 4, center = FALSE)$a, fit(x, center = FALSE)$a,
               tolerance = 1e-6)
})

test_that("spearman uses the location and scale of the data it is applied to", {
  skip_if_no_torch()
  x <- make_data()
  fit <- function(x, ...) fit_quick(x, transformation = "spearman", epochs = 30, ...)
  # scaled but not centered: invariant to scaling of the data
  expect_equal(fit(x * 4, center = FALSE)$a, fit(x, center = FALSE)$a,
               tolerance = 1e-6)
  # without standardization, groups play no role for the initial estimate
  expect_equal(fit(x, groups = rep(1:2, 20), center = FALSE, scale = FALSE)$a,
               fit(x, center = FALSE, scale = FALSE)$a, tolerance = 1e-6)
})

test_that("with groups, location and scale of every group are returned", {
  skip_if_no_torch()
  x <- make_data()
  groups <- rep(c("b", "a"), 20)
  res <- fit_quick(x, groups = groups, epochs = 30)
  expect_equal(dim(res$robcenter), c(2, 6))
  expect_equal(dim(res$robscale), c(2, 6))
  expect_equal(rownames(res$robscale), c("b", "a"))
  z <- x
  for (g in unique(groups)) {
    rows <- groups == g
    z[rows, ] <- scale(x[rows, ], center = res$robcenter[g, ],
                       scale = res$robscale[g, ])
  }
  expect_equal(res$phi, z %*% res$a, ignore_attr = TRUE)
  # nothing is estimated without centering / scaling
  res_raw <- fit_quick(x, groups = groups, center = FALSE, scale = FALSE,
                       epochs = 30)
  expect_true(all(res_raw$robcenter == 0))
  expect_true(all(res_raw$robscale == 1))
})
