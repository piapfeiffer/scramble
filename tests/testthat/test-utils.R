test_that("principal_angles is zero for identical and pi/2 for orthogonal subspaces", {
  a <- diag(4)[, 1:2]
  expect_equal(as.numeric(principal_angles(a, a)$angles), c(0, 0))
  expect_equal(as.numeric(principal_angles(a, a %*% matrix(c(0, 1, 1, 0), 2))$angles),
               c(0, 0))
  b <- diag(4)[, 3:4]
  expect_equal(as.numeric(principal_angles(a, b)$angles), c(pi / 2, pi / 2))
})

test_that("principal_angles recovers a known angle", {
  a <- matrix(c(1, 0, 0), 3)
  b <- matrix(c(1, 1, 0) / sqrt(2), 3)
  expect_equal(as.numeric(principal_angles(a, b)$angles), pi / 4)
})

test_that("TPR and TNR compare sparsity patterns", {
  A_true <- matrix(c(1, 1, 0, 0), 2)
  A <- matrix(c(2, 0, 0, 3), 2)
  expect_equal(TPR(A_true, A), 0.5)
  expect_equal(TNR(A_true, A), 0.5)
  expect_equal(TPR(A_true, A_true), 1)
  expect_equal(TNR(A_true, A_true), 1)
  expect_visible(TPR(A_true, A))
})

test_that("f1_score counts true and false positives", {
  expect_equal(f1_score(1:4, 1:4), 1)
  expect_equal(f1_score(1:4, 5:8), 0)
  # TP = 2, FP = 1, FN = 2
  expect_equal(f1_score(1:4, c(1, 2, 9)), 2 / (2 + 0.5 * 3))
})

test_that("penalty_fun_unc is the elastic net penalty", {
  skip_if_no_torch()
  z <- c(0.5, -0.3, 0, 0.8)
  zt <- torch::torch_tensor(z, dtype = torch::torch_double())
  expected <- function(alpha) 2 * (alpha * sum(abs(z)) + (1 - alpha) * sum(z^2))
  for (alpha in c(0, 0.5, 1)) {
    expect_equal(as.numeric(penalty_fun_unc(zt, 2, alpha)), expected(alpha),
                 tolerance = 1e-6)
  }
})

test_that("robust losses return non-negative scalars and are zero at zero", {
  skip_if_no_torch()
  m <- torch::torch_tensor(matrix(rnorm(30), 10, 3), dtype = torch::torch_double())
  for (loss in list(loss_huber, loss_tukey, loss_lts)) {
    val <- loss(m)
    expect_equal(val$numel(), 1)
    expect_gte(as.numeric(val), 0)
  }
  zero <- torch::torch_zeros(5, 2, dtype = torch::torch_double()) + 1e-12
  expect_equal(as.numeric(loss_huber(zero)), 0, tolerance = 1e-10)
  expect_equal(as.numeric(loss_lts(zero)), 0, tolerance = 1e-10)
})

test_that("Huber loss downweights an outlier compared to squared loss", {
  skip_if_no_torch()
  r <- c(rnorm(19, sd = 0.5), 50)
  m <- torch::torch_tensor(matrix(r, 20, 1), dtype = torch::torch_double())
  expect_lt(as.numeric(loss_huber(m)), mean(r^2) / 10)
})
