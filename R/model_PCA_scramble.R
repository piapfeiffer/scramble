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
    self$TV <- apply(X, 2, robustbase::Qn)
    # converting X to a tensor is costly for large data; do it only once.
    # The residuals have always been computed as `X - B` with the R matrix X,
    # which torch converts with its default dtype (float32); X_resid keeps
    # exactly that behaviour so that results are unchanged.
    self$X_tensor <- torch::torch_tensor(X, dtype = torch::torch_double())
    self$X_resid <- torch::torch_tensor(X)

  }, forward = function(X) {
    # X must be the matrix the module was initialized with; its tensor
    # versions are cached in initialize()
    X_t <- self$X_tensor
    X_r <- self$X_resid

    m1 <- torch::torch_tensor(rep(0.0, ncol(self$V)),
                              dtype = torch::torch_double())
    for (j in 1:ncol(self$V)){
      m1[j] <- penalty_fun_unc(self$V[,j], self$pen_x * self$TV[j], self$alpha_x)
    }

    G <- torch::torch_cat(c(torch::torch_unsqueeze(m1, 2)))

    B <-  torch::torch_matmul(X_t,
                              torch::torch_matmul(self$V,
                                                  torch::torch_transpose(self$V, 1, -1)))
    if (self$loss_type == "L2") {
      loss <- torch::torch_mean(torch::torch_square(X_r - B)) + torch::torch_sum(G)

    } else if (self$loss_type == "Huber") {
      loss <- loss_huber(X_r - B) + torch::torch_sum(G)

    } else if (self$loss_type == "Tukey") {
      loss <- loss_tukey(X_r - B) + torch::torch_sum(G)

    } else if (self$loss_type == "LTS") {
      if (is.na(self$groups[1])){
        loss <- loss_lts(X_r - B, self$h) + torch::torch_sum(G)
      }
      else{
        loss <- torch::torch_sum(G)
        for (group in unique(self$groups)){
          XG <- X[which(self$groups == group),]
          BG <-  torch::torch_matmul(torch::torch_tensor(XG, dtype = torch::torch_double()),
                                    torch::torch_matmul(self$V,
                                                        torch::torch_transpose(self$V, 1, -1)))
          loss <- loss + loss_lts(XG - BG, self$h)
        }
      }

    }else {
      warning(paste("Loss type:",self$loss_type, "not implemented, standard L2 loss is used instead!"))
      loss <- torch::torch_square(torch::torch_norm(X_r - B)) + torch::torch_sum(G)
    }

     loss
  }, evaluate = function(X, V) {
    res <- apply(X %*% V, 2, robustbase::Qn)
    return(res)
  }
)
