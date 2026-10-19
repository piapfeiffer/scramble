# the objective of model_PCA_SVD, computed in R in double precision
objective_r <- function(x, v, pen_x, alpha_x, loss, groups = NA) {
  resid <- x - x %*% v %*% t(v)
  penalty <- pen_x * apply(x, 2, robustbase::Qn)[seq_len(ncol(v))] *
    (alpha_x * colSums(tanh(1000 * v) * v) + (1 - alpha_x) * colSums(v^2))
  if (is.na(groups[1])) {
    value <- loss(resid)
  } else {
    value <- sum(vapply(split(seq_len(nrow(x)), groups),
                        function(rows) loss(resid[rows, , drop = FALSE]),
                        numeric(1)))
  }
  value + sum(penalty)
}

lts_r <- function(resid, h = 0.75) {
  mean(apply(resid^2, 2, function(r) {
    threshold <- sort(r)[round((length(r) - 1) * h) + 1]
    ifelse(r <= threshold, r, 0)
  }))
}

test_that("the objective is computed in double precision", {
  skip_if_no_torch()
  x <- scale(make_data())
  v <- svd(x)$v[, 1:2]
  groups <- rep(1:3, length.out = nrow(x))
  model <- function(loss_type, groups = NA) {
    model_PCA_SVD(x, groups, list(v = v), nrow(x), ncol(x), pen_x = 0.1,
                  alpha_x = 0.5, loss_type = loss_type, h = 0.75)
  }
  expect_equal(as.numeric(model("L2")(x)),
               objective_r(x, v, 0.1, 0.5, function(r) mean(r^2)),
               tolerance = 1e-12)
  expect_equal(as.numeric(model("LTS")(x)),
               objective_r(x, v, 0.1, 0.5, lts_r), tolerance = 1e-12)
  expect_equal(as.numeric(model("LTS", groups)(x)),
               objective_r(x, v, 0.1, 0.5, lts_r, groups), tolerance = 1e-12)
})
