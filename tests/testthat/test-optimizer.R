test_that("the manifold optimizer keeps the loadings orthonormal", {
  skip_if_no_torch()
  set.seed(2)
  v0 <- qr.Q(qr(matrix(rnorm(12), 6, 2)))
  V <- torch::torch_tensor(v0, dtype = torch::torch_double(), requires_grad = TRUE)
  target <- torch::torch_tensor(matrix(rnorm(12), 6, 2), dtype = torch::torch_double())
  # fixed step size without line search: lr must be small enough not to overshoot
  opt <- optim_sgd_manifold(V, lr = 0.01, manifold = TRUE)
  losses <- numeric(20)
  for (i in 1:20) {
    opt$zero_grad()
    loss <- torch::torch_sum((V - target)^2)
    loss$backward()
    opt$step()
    losses[i] <- as.numeric(loss)
  }
  expect_true(is_orthonormal(as.matrix(V$detach())))
  expect_false(is.unsorted(rev(losses)))
})

test_that("the optimizer does Euclidean steps when manifold = FALSE", {
  skip_if_no_torch()
  V <- torch::torch_tensor(c(1, 2), dtype = torch::torch_double(), requires_grad = TRUE)
  opt <- optim_sgd_manifold(V, lr = 0.5, manifold = FALSE)
  loss <- torch::torch_sum(V^2)
  loss$backward()
  opt$step()
  expect_equal(as.numeric(V$detach()), c(0, 0))
})
