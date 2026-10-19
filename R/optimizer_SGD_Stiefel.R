# Riemannian SGD on the Stiefel manifold (Cayley transform retraction).
#
# The structure of this optimizer (initialize/step, the loop over parameter
# groups and parameters, and the plain Euclidean SGD step for
# manifold = FALSE) is adapted from optim_sgd() and the optimizer template in
# the 'torch' package, Copyright (c) 2020 Daniel Falbel and RStudio, licensed
# under the MIT license. See inst/COPYRIGHTS.
optim_sgd_manifold <- torch::optimizer(
  initialize = function(params, lr, manifold, threshold = NA) {
    defaults <- list(
      lr = lr,
      manifold = manifold,
      threshold = threshold
    )
    super$initialize(params, defaults)
  },
  step = function() {

    torch::with_no_grad({

        for (g in seq_along(self$param_groups)) {
          group <- self$param_groups[[g]]
          if (isTRUE(group$manifold)){

            for (p in seq_along(group$params)) {
              param <- group$params[[p]]

              if (is.null(param$grad) || torch::is_undefined_tensor(param$grad)) {
                next
              }

              manifold_grad <- -group$lr * (param$grad - torch::torch_matmul(param,
                                                                torch::torch_matmul(
                                                                  torch::torch_transpose(param$grad, 1, -1),
                                                                  param)
                                                                ))

              # Cayley transform (I - A/2)^{-1} (I + A/2) param with the
              # rank-2k matrix A = W param^T - param W^T, solved via the
              # Woodbury identity in a (2k x 2k) system instead of (p x p)
              U <- torch::torch_cat(list(manifold_grad, param), dim = 2)
              Ct <- torch::torch_transpose(
                torch::torch_cat(list(param, -manifold_grad), dim = 2), 1, -1)
              rhs <- param + 1/2 * torch::torch_matmul(U, torch::torch_matmul(Ct, param))
              small <- torch::torch_eye(U$size(2), dtype = torch::torch_double()) -
                1/2 * torch::torch_matmul(Ct, U)
              new_param <- rhs + 1/2 * torch::torch_matmul(
                U, torch::linalg_solve(small, torch::torch_matmul(Ct, rhs)))

              if (!is.na(group$threshold)){
                new_param <- soft(new_param, group$threshold)
                res_svd <- torch::linalg_svd(new_param, full_matrices=FALSE)
                proj_param <- torch::torch_matmul(res_svd[[1]], res_svd[[3]])
                param$copy_(proj_param)
              } else {
                param$copy_(new_param)
              }
            }
          }
       else {
          for (p in seq_along(group$params)) {
            param <- group$params[[p]]

            if (is.null(param$grad) || torch::is_undefined_tensor(param$grad)) {
              next
            }

            param$add_(param$grad, alpha = -group$lr)
          }
        }
      }

    })
  }
)
