#' Torch module for robust and sparse PCA
#'
#' The loadings `V` are the only parameter. The forward pass returns the
#' (robust) reconstruction loss of `X` from its projection onto `V` plus the
#' elastic net penalty of each column of `V`.
#'
#' @param X Standardized (n x p) data matrix.
#' @param groups Vector of group labels or `NA`; only used by the LTS loss.
#' @param rsvd List with element `v`, the (p x k) initial loadings.
#' @param n,p Number of observations and variables.
#' @param pen_x Penalty parameter.
#' @param alpha_x Elastic net mixing parameter.
#' @param loss_type One of "L2", "Huber", "Tukey", "LTS".
#' @param h Proportion of cells kept per variable for the LTS loss.
#' @param warm_v Optional warm start for the loadings, used instead of `rsvd$v`.
#' @param warm_m Currently unused.
#' @noRd
model_PCA_SVD <- torch::nn_module(
  initialize = function(X,
                        groups,
                        rsvd,
                        n,
                        p,
                        pen_x,
                        alpha_x,
                        loss_type,
                        h,
                        warm_v = NA,
                        warm_m = NA) {
    self$loss_type <- loss_type
    self$p <- p
    self$n <- n
    self$h <- h
    self$pen_x <- pen_x
    self$groups <- groups

    self$alpha_x <- alpha_x

    if (is.na(warm_v[1])) {
        v <- rsvd$v
        self$V <- torch::nn_parameter(torch::torch_tensor(v,
                                                          dtype = torch::torch_double()), requires_grad = TRUE)

      } else {
        v <- warm_v
        self$V <- torch::nn_parameter(torch::torch_tensor(v,
                                                          dtype = torch::torch_double()), requires_grad = TRUE)
      }
    # The penalty of component j is weighted by the Qn of *variable* j. With
    # the default robust scaling (scale = TRUE, no groups) all these weights
    # are 1. Kept as is deliberately, to reproduce published results.
    self$TV <- apply(X, 2, robustbase::Qn)
    # converting X to a tensor is costly for large data; do it only once.
    # Before version 0.2.0 the residuals were computed from a float32 copy of
    # X (torch's default dtype); now everything is in double precision.
    self$X_tensor <- torch::torch_tensor(X, dtype = torch::torch_double())
    if (!is.na(groups[1])) {
      self$group_rows <- lapply(unique(groups), function(g) which(groups == g))
    }

  }, forward = function(X) {
    # X must be the matrix the module was initialized with; its tensor
    # versions are cached in initialize()
    X_t <- self$X_tensor

    # elastic net penalty of all columns at once
    V <- self$V
    l1 <- torch::torch_sum(torch::torch_tanh(1000 * V) * V, dim = 1)
    l2 <- torch::torch_sum(V^2, dim = 1)
    weights <- torch::torch_tensor(self$pen_x * self$TV[seq_len(ncol(V))],
                                   dtype = torch::torch_double())
    G <- (self$alpha_x * l1 + (1 - self$alpha_x) * l2) * weights

    B <- torch::torch_matmul(torch::torch_matmul(X_t, V),
                             torch::torch_transpose(V, 1, -1))
    resid <- X_t - B
    if (self$loss_type == "L2") {
      loss <- torch::torch_mean(torch::torch_square(resid)) + torch::torch_sum(G)

    } else if (self$loss_type == "Huber") {
      loss <- loss_huber(resid) + torch::torch_sum(G)

    } else if (self$loss_type == "Tukey") {
      loss <- loss_tukey(resid) + torch::torch_sum(G)

    } else if (self$loss_type == "LTS") {
      if (is.na(self$groups[1])) {
        loss <- loss_lts(resid, self$h) + torch::torch_sum(G)
      } else {
        # trimming per group: sum of the LTS losses of the groups
        loss <- torch::torch_sum(G)
        for (rows in self$group_rows) {
          loss <- loss + loss_lts(resid[rows, , drop = FALSE], self$h)
        }
      }

    } else {
      warning(paste("Loss type:",self$loss_type, "not implemented, standard L2 loss is used instead!"))
      loss <- torch::torch_square(torch::torch_norm(resid)) + torch::torch_sum(G)
    }

     loss
  }, evaluate = function(X, V) {
    res <- apply(X %*% V, 2, robustbase::Qn)
    return(res)
  }
)
