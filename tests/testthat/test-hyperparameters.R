test_that("the penalty is tuned by Bayesian optimization", {
  skip_if_no_torch()
  skip_on_cran()
  set.seed(42)
  res <- suppressWarnings(pcaSCRAMBLE(make_data(), k = 2, alpha_x = 0.5,
                                      epochs = 2000, tol = 1e-4))
  expect_gte(res$pen_x, 1e-3)
  expect_lte(res$pen_x, 10)
  expect_s3_class(res$summary, "data.frame")
  expect_true(all(c("pen_x", "Score", "Pred", "Sparse_x") %in% names(res$summary)))
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(res))
})

test_that("a failed optimization falls back to the upper bound with a warning", {
  skip_if_no_torch()
  # with so few epochs no fit converges, all scores are 0 and the GP fails
  set.seed(1)
  out <- collect_warnings(pcaSCRAMBLE(make_data(), k = 2, alpha_x = 0.5,
                                      epochs = 5, bounds = c(0.01, 3)))
  expect_true(any(grepl("Hyperparameter optimization failed", out$warnings)))
  expect_equal(out$value$pen_x, 3)
  expect_false(anyNA(out$value$a))
})
