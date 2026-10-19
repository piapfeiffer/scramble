#' Train the robust and sparse PCA model
#'
#' Low-level routine used by [pcaSCRAMBLE()]: runs Riemannian gradient descent
#' on the Stiefel manifold for a fixed penalty parameter, starting from the
#' loadings in `rsvd$v`, and thresholds small loadings afterwards. Most users
#' should call [pcaSCRAMBLE()] instead, which also standardizes the data and
#' tunes the penalty.
#'
#' @param X An (n x p) numeric data matrix, already standardized.
#' @param groups Vector of length n with the group of each observation, or
#'  `NA` for no groups. Only used by the `"LTS"` loss.
#' @param rsvd A list with element `v`, a (p x k) matrix of initial loadings
#'  (e.g. the result of [svd()] restricted to k columns).
#' @param n Number of observations.
#' @param p Number of variables.
#' @param pen_x Penalty parameter (non-negative scalar).
#' @param alpha_x Elastic net mixing parameter between 0 (ridge) and 1 (lasso).
#' @param loss_type One of `"L2"`, `"Huber"`, `"Tukey"` or `"LTS"`.
#' @param h Proportion of cells kept per variable for the `"LTS"` loss.
#' @param tol Convergence tolerance on the change of the objective.
#' @param lr Learning rate.
#' @param epochs Maximum number of epochs.
#' @param lr_decay Multiplicative learning rate decay per epoch.
#' @param show_warnings Logical: warn if the algorithm did not converge? If
#'  `FALSE`, `best_measure` is set to 0 instead.
#' @param init_ortho Currently unused.
#' @returns A list with elements `best_measure` (robust variances of the sparse
#'  components), `best_loss`, `best_v` (sparse loadings), `best_v_ortho`
#'  (orthogonal loadings before thresholding), `plot_loss` (loss per epoch)
#'  and `best_measure_ortho`.
#' @export
#' @examples
#' if (torch::torch_is_installed()) {
#'   set.seed(1)
#'   x <- scale(matrix(rnorm(40 * 5), 40, 5))
#'   rsvd <- svd(x)
#'   rsvd$v <- rsvd$v[, 1:2]
#'   res <- train_PCA_SVD(x, groups = NA, rsvd = rsvd, n = 40, p = 5,
#'                        pen_x = 0.1, alpha_x = 0.5, loss_type = "L2",
#'                        h = NA, tol = 1e-5, lr = 1e-3, epochs = 100,
#'                        lr_decay = 1)
#'   res$best_v
#' }
train_PCA_SVD <- function(X,
                          groups,
                          rsvd, n, p,
                          pen_x,
                          alpha_x,
                          loss_type,
                          h,
                          tol, lr, epochs,
                          lr_decay,
                          show_warnings = FALSE,
                          init_ortho = TRUE) {

  plot_loss <- rep(NA, epochs)

  model <- model_PCA_SVD(X,
                         groups,
                         rsvd, n, p,
                         pen_x = pen_x,
                         alpha_x = alpha_x,
                         loss_type = loss_type,
                         h = h,
                         warm_v = NA,
                         warm_m = NA)

  best_v <- as.matrix(model$parameters$V$detach())

  best_loss <- 1e16
  improved <- FALSE
  change <- 1e16
  count_ep <- 0

  old_v <- best_v
  old_loss <- best_loss

  optimizer2 <- optim_sgd_manifold(model$parameters$V, lr = lr, manifold = TRUE)
  lmbda <- function(epoch) lr_decay
  scheduler2 <- torch::lr_multiplicative(optimizer2, lr_lambda = lmbda)

  model$zero_grad()

  # relative changes of the loadings in the 10 most recent epochs (newest
  # first), used for the thresholding after convergence
  recent_diff_v <- numeric(0)

  while (as.numeric(change) > tol && count_ep < epochs) {
    optimizer2$zero_grad()

    loss <- model(X)
    plot_loss[count_ep + 1] <- as.numeric(loss)

    if (as.numeric(loss) < best_loss) {
      best_loss <- as.numeric(loss)
      best_v <- as.matrix(model$parameters$V$detach())
      improved <- TRUE
    }

    loss$backward()

    optimizer2$step()
    scheduler2$step()

    new_v <- as.matrix(model$parameters$V$detach())
    if (norm(old_v, "1") != 0) {
      diff_v <- norm(new_v - old_v, "F") / norm(old_v, "F")
    } else {
      diff_v <- 0
    }
    recent_diff_v <- head(c(diff_v, recent_diff_v), 10)

    old_v <- new_v

    count_ep <- count_ep + 1
    change <- abs(as.numeric(old_loss - loss))
    old_loss <- loss

  } # end inner loop

  # robust variances of the best loadings (computed once, not every epoch)
  best_measure <- if (improved) model$evaluate(X, best_v) else 0

  norm_vec <- function(x) sqrt(sum(x^2))

  # set loadings to zero that are smaller than the recent fluctuation of the
  # iterates, then renormalize the columns
  # (the sd needs two epochs: after a single epoch only the mean is used,
  # without any epoch nothing is thresholded)
  threshold_v <- switch(min(length(recent_diff_v), 2) + 1,
                        0,
                        recent_diff_v,
                        mean(recent_diff_v) + 2 * sd(recent_diff_v))
  best_v_ortho <- best_v
  best_measure_ortho <- best_measure
  best_v_thresh <- ifelse(abs(best_v) > threshold_v, best_v, 0)
  if (norm(as.matrix(best_v_thresh), "1") > 0) {
    best_v <- best_v_thresh
    col_norms <- apply(best_v, 2, norm_vec)
    best_v <- scale(best_v,
                    center = FALSE,
                    scale = ifelse(col_norms == 0, 1, col_norms))
    best_measure <- model$evaluate(X, best_v)
  }

  if (as.numeric(change) > tol) {
    if (isTRUE(show_warnings)) {
      warning(paste("Loop did not converge. Try more iterations. Change:",
                    as.numeric(change)))
    } else {
      best_measure <- 0
    }
  }

  list(best_measure = best_measure,
       best_loss = best_loss,
       best_v = best_v,
       best_v_ortho = best_v_ortho,
       plot_loss = plot_loss,
       best_measure_ortho = best_measure_ortho)
}
