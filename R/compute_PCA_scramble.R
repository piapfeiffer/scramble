#' Cellwise robust and sparse PCA
#'
#' Computes sparse and cellwise robust principal components by Riemannian
#' gradient descent on the Stiefel manifold. The data matrix is decomposed
#' directly instead of forming the full covariance matrix, which keeps the
#' memory footprint small even for many variables.
#'
#' @param data_x An (n x p) numeric data matrix (or data frame).
#' @param groups Optional vector of length n containing the group of each
#'  observation. If `NA` (default), no group structure is used. With groups,
#'  centering and scaling are done per group.
#' @param transformation Data transformation used for the initial estimate,
#'  one of `"identity"` (default, corresponds to the sample covariance),
#'  `"spearman"` (rank-based) or `"wrapping"`.
#' @param loss_type Loss function, one of `"L2"` (default, least squares),
#'  `"Huber"`, `"Tukey"` or `"LTS"` (least trimmed squares).
#' @param param Parameter of the `"LTS"` loss: the proportion of cells that is
#'  kept per variable. Try 0.75 (used if `NA`) and increase or decrease for
#'  tuning. Ignored for the other loss types.
#' @param center Logical: whether the data should be robustly centered (using
#'  the median) before proceeding. Default is `TRUE`.
#' @param scale Logical: whether the data should be robustly scaled (using the
#'  Qn estimator, or the MAD if `groups` are given) before proceeding.
#'  Default is `TRUE`.
#' @param alpha_x Elastic net mixing parameter between 0 (ridge) and 1 (lasso).
#'  Must be supplied.
#' @param k Number of principal components to compute. Defaults to
#'  `min(n, p)`; a smaller value speeds up computation considerably.
#' @param tol Convergence tolerance of the gradient descent.
#' @param lr Learning rate of the gradient descent.
#' @param epochs Maximum number of gradient descent iterations.
#' @param lr_decay Multiplicative learning rate decay per epoch.
#' @param criterion Criterion for selecting the penalty parameter: `"TPO"`
#'  (default, tradeoff product optimization), `"BIC"`, or any other value to
#'  maximize the explained variance only.
#' @param bounds Numeric vector of length 2, the search range of the penalty
#'  parameter used in the Bayesian optimization.
#' @param penalties Optional list with element `pen_x`. If given, no
#'  hyperparameter optimization is done and `pen_x` is used as penalty.
#'
#' @returns An object of class `"PCA_result"`, a list with elements
#' \describe{
#'   \item{a}{(p x k) matrix of sparse loadings.}
#'   \item{a_ortho}{(p x k) matrix of orthogonal loadings before thresholding.}
#'   \item{u, s}{Currently always `NULL`; kept for backward compatibility.}
#'   \item{loss}{Value of the objective in each epoch (`NA` after convergence).}
#'   \item{measure}{Robust variances (Qn) of the sparse components.}
#'   \item{measure_ortho}{Robust variances (Qn) of the orthogonal components.}
#'   \item{explained_var}{Cumulative proportion of explained variance.}
#'   \item{phi}{(n x k) matrix of scores, the standardized data times `a`.}
#'   \item{pen_x}{The penalty parameter used.}
#'   \item{summary}{Summary of the hyperparameter optimization, if done.}
#'   \item{wt}{Currently always `NA`; kept for backward compatibility.}
#'   \item{robscale, robcenter}{Robust scale and location estimates of the
#'     variables. With `groups`, a matrix with one row per group (row names
#'     are the group labels) containing the estimates used for that group;
#'     0 (location) and 1 (scale) if `center` or `scale` is `FALSE`.}
#'   \item{alpha_x}{The elastic net parameter used.}
#' }
#' @export
#' @examples
#' if (torch::torch_is_installed()) {
#'   set.seed(1)
#'   x <- matrix(rnorm(50 * 6), 50, 6)
#'   x[, 2] <- x[, 1] + rnorm(50, sd = 0.1)
#'   # fixed penalty: no hyperparameter optimization
#'   res <- pcaSCRAMBLE(x, alpha_x = 0.5, k = 2,
#'                      penalties = list(pen_x = 0.1), epochs = 200)
#'   res$a
#'   plot(res)
#' }
#' \donttest{
#' if (torch::torch_is_installed()) {
#'   # robust version with tuned penalty (fewer epochs to keep it short)
#'   res_rob <- pcaSCRAMBLE(x, transformation = "wrapping", loss_type = "Huber",
#'                          alpha_x = 0.5, k = 2, epochs = 500, tol = 1e-4)
#' }
#' }
pcaSCRAMBLE <- function(data_x,
                        groups = NA,
                        transformation = "identity",
                        loss_type = "L2",
                        param = NA,
                        center = TRUE,
                        scale = TRUE,
                        alpha_x = NA,
                        k = NA,
                        tol = 1e-5,
                        lr = 1e-3,
                        epochs = 2000,
                        lr_decay = 1,
                        criterion = "TPO",
                        bounds = c(1e-3, 10),
                        penalties = NA) {
  # note: the argument `scale` shadows base::scale(), which is still found
  # below because R skips non-function objects when looking up a call.

  data_x <- as.matrix(data_x)
  if (!is.numeric(data_x)) {
    stop("`data_x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  p <- ncol(data_x)
  n <- nrow(data_x)

  if (!transformation %in% c("identity", "spearman", "wrapping")) {
    stop("`transformation` must be one of \"identity\", \"spearman\", \"wrapping\".",
         call. = FALSE)
  }
  if (!loss_type %in% c("L2", "Huber", "Tukey", "LTS")) {
    stop("`loss_type` must be one of \"L2\", \"Huber\", \"Tukey\", \"LTS\".",
         call. = FALSE)
  }
  if (!is.numeric(alpha_x) || length(alpha_x) != 1 || is.na(alpha_x) ||
      alpha_x < 0 || alpha_x > 1) {
    stop("`alpha_x` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.na(groups[1]) && length(groups) != n) {
    stop("`groups` must have one entry per row of `data_x`.", call. = FALSE)
  }
  if (loss_type == "LTS") {
    if (is.na(param)) {
      param <- 0.75
    }
    if (!is.numeric(param) || param <= 0 || param > 1) {
      stop("`param` must be in (0, 1] for loss_type = \"LTS\".", call. = FALSE)
    }
  }
  if (!is.numeric(epochs) || length(epochs) != 1 || is.na(epochs) || epochs < 1) {
    stop("`epochs` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(bounds) || length(bounds) != 2 || any(bounds < 0)) {
    stop("`bounds` must be a numeric vector of two non-negative values.",
         call. = FALSE)
  }

  if (is.na(k)) {
    # per default, compute best possible rank k approximation
    # in case of really big matrices, it can speed up computation to provide a
    # smaller number for an initial result
    k <- min(p, n)
  }
  if (k < 1) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  if (k > min(p, n)) {
    warning("k has to be <= min(p, n) ! k = min(p, n) is used")
    k <- min(p, n)
  }

  wt <- rep(NA, n)

  robcenter <- apply(data_x, 2, median, na.rm = TRUE)
  robscale <- apply(data_x, 2, robustbase::Qn, na.rm = TRUE)

  # first, the data is standardized
  if (is.na(groups[1])) { # no groups
    if (isTRUE(center)) {
      data_x <- scale(data_x,
                      center = robcenter,
                      scale = FALSE)
    }
    if (isTRUE(scale)) {
      data_x <- scale(data_x,
                      center = FALSE,
                      scale = ifelse(robscale == 0, 1, robscale))
    }

  } else { # groups
    # location and scale used for each group, one row per group (before
    # version 0.2.0 only those of the last group were returned)
    group_labels <- unique(groups)
    robcenter <- matrix(0, length(group_labels), p,
                        dimnames = list(group_labels, colnames(data_x)))
    robscale <- matrix(1, length(group_labels), p,
                       dimnames = list(group_labels, colnames(data_x)))
    for (i in seq_along(group_labels)) {
      rows <- which(groups == group_labels[i])
      data_g <- data_x[rows, , drop = FALSE]

      if (isTRUE(center)) {
        robcenter[i, ] <- apply(data_g, 2, median, na.rm = TRUE)
        data_g <- scale(data_g,
                        center = robcenter[i, ],
                        scale = FALSE)
        data_x[rows, ] <- data_g
      }

      if (isTRUE(scale)) {
        robscale[i, ] <- apply(data_g, 2, stats::mad, na.rm = TRUE)
        data_g <- scale(data_g,
                        center = FALSE,
                        scale = ifelse(robscale[i, ] == 0, 1, robscale[i, ]))
        data_x[rows, ] <- data_g
      }
    }
  }

  # robust location and scale of the data entering the transformation: data
  # standardized above have location 0 and scale 1
  loc_x <- if (isTRUE(center)) rep(0, p) else apply(data_x, 2, median, na.rm = TRUE)
  scale_x <- if (isTRUE(scale)) rep(1, p) else apply(data_x, 2, robustbase::Qn, na.rm = TRUE)

  # robust initial estimate of the loadings
  if (transformation == "identity") {
    rsvd <- svd(data_x)

  } else if (transformation == "spearman") {
    t_x <- apply(data_x, 2, rank, ties.method = "average")
    t_x <- scale(t_x,
                 center = apply(t_x, 2, mean),
                 scale = apply(t_x, 2, sd))
    # give the ranks the location and scale of the data; before version
    # 0.2.0 those of the raw data (or, with groups, of the last group) were used
    if (!isTRUE(scale)) {
      t_x <- sweep(t_x, MARGIN = 2, scale_x, `*`)
    }
    if (!isTRUE(center)) {
      t_x <- sweep(t_x, MARGIN = 2, loc_x, `+`)
    }
    rsvd <- svd(t_x)

  } else if (transformation == "wrapping") {
    # before version 0.2.0, the location and scale of the raw data were used
    # here, which standardized the data a second time
    t_x <- cellWise::wrap(data_x,
                          locX = loc_x,
                          scaleX = scale_x)$Xw
    rsvd <- svd(t_x)
  }

  # order the initial components by their robust variance and keep k of them
  res <- apply(data_x %*% rsvd$v, 2, robustbase::Qn)
  order_res <- order(res, decreasing = TRUE)
  rsvd$v <- rsvd$v[, order_res[1:k], drop = FALSE]

  SUMMARY <- list()
  if (any(is.na(penalties))) {
    res_param <- bayesian_optimization_PCA_SVD(data_x, groups,
                                               rsvd, n, p, k,
                                               alpha_x = alpha_x,
                                               loss_type = loss_type,
                                               h = param,
                                               tol = tol,
                                               lr = lr,
                                               epochs = epochs,
                                               lr_decay = lr_decay,
                                               bounds_input = bounds,
                                               criterion = criterion)
    PEN_X <- as.numeric(res_param$best_params)
    SUMMARY <- res_param$summary
  } else {
    PEN_X <- penalties$pen_x
  }

  res <- train_PCA_SVD(data_x,
                       groups,
                       rsvd,
                       n,
                       p,
                       pen_x = PEN_X,
                       alpha_x = alpha_x,
                       loss_type = loss_type,
                       h = param,
                       tol = tol,
                       lr = lr,
                       epochs = epochs,
                       lr_decay = lr_decay,
                       show_warnings = TRUE)

  CORR <- res$best_measure
  CORR_ortho <- res$best_measure_ortho
  order_corr <- order(CORR, decreasing = TRUE)
  order_corr_ortho <- order(CORR_ortho, decreasing = TRUE)

  A <- res$best_v[, order_corr]
  A_ortho <- res$best_v_ortho[, order_corr_ortho]
  PHI <- as.matrix(data_x) %*% as.matrix(A)
  measure_ortho <- as.numeric(CORR_ortho)[order_corr_ortho]

  structure(list(
    a = A,
    a_ortho = A_ortho,
    u = NULL, # not computed, kept for backward compatibility
    s = NULL, # not computed, kept for backward compatibility
    loss = res$plot_loss,
    measure = as.numeric(CORR)[order_corr],
    measure_ortho = measure_ortho,
    explained_var = cumsum(measure_ortho) / sum(measure_ortho),
    phi = PHI,
    pen_x = PEN_X,
    summary = SUMMARY,
    wt = wt,
    robscale = robscale,
    robcenter = robcenter,
    alpha_x = alpha_x
  ), class = "PCA_result")
}

#' Plot a PCA result
#'
#' Plots the explained variance, the robust variances, the scores of the first
#' two components, the loss trace and, if the penalty was tuned, the
#' hyperparameter optimization.
#'
#' @param x An object of class `"PCA_result"`, as returned by [pcaSCRAMBLE()].
#' @param ... Currently unused.
#' @returns `x`, invisibly.
#' @export
#' @examples
#' if (torch::torch_is_installed()) {
#'   set.seed(1)
#'   x <- matrix(rnorm(50 * 6), 50, 6)
#'   res <- pcaSCRAMBLE(x, alpha_x = 0.5, k = 2,
#'                      penalties = list(pen_x = 0.1), epochs = 200)
#'   plot(res)
#' }
plot.PCA_result <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(mfrow = c(1, 3))
  ev <- round(x$explained_var, 2)
  graphics::barplot(ev, xlab = "# PC", ylab = "Proportion of explained variance")
  plot(x$measure_ortho, xlab = "# PC", ylab = "Variances", type = "b", lty = 2, pch = 1)
  if (NCOL(x$phi) >= 2) {
    plot(x$phi[, 1], x$phi[, 2], xlab = "PC 1", ylab = "PC 2", type = "p", pch = 1)
  } else {
    plot(x$phi[, 1], xlab = "Observation", ylab = "PC 1", type = "p", pch = 1)
  }

  graphics::par(mfrow = c(1, 3))
  plot(x$loss, type = "l", xlab = "Epoch", ylab = "Loss")
  if (is.data.frame(x$summary) && nrow(x$summary) > 0) {
    ord <- order(x$summary$pen_x)
    plot(x$summary$pen_x[ord], x$summary$Pred[ord], type = "l",
         xlab = "pen_x", ylab = "Total variance")
    graphics::abline(v = x$pen_x[1], lty = 2, col = "red")
    plot(x$summary$pen_x[ord], x$summary$Sparse_x[ord], type = "l",
         xlab = "pen_x", ylab = "Sparsity")
    graphics::abline(v = x$pen_x[1], lty = 2, col = "red")
  }
  invisible(x)
}
