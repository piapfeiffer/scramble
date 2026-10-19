#' Principal angles between subspaces
#'
#' Calculates the principal angles between the column spaces of `a` and `b`,
#' e.g. to compare estimated and true loadings.
#'
#' @param a A matrix.
#' @param b A matrix with the same number of rows as `a`.
#' @returns A list with element `angles`, a column matrix of the principal
#'  angles (in radians).
#' @export
#' @examples
#' a <- diag(3)[, 1:2]
#' b <- cbind(c(1, 1, 0) / sqrt(2), c(0, 0, 1))
#' principal_angles(a, b)
principal_angles <- function(a, b) {
  angles <- matrix(0, ncol = ncol(a), nrow = 1)
  qa <- qr.Q(qr(a))
  qb <- qr.Q(qr(b))
  C <- svd(t(qa) %*% qb)$d
  rkA <- qr(a)$rank
  rkB <- qr(b)$rank
  if (rkA <= rkB) {
    B <- qb - qa %*% (t(qa) %*% qb)
  } else {
    B <- qa - qb %*% (t(qb) %*% qa)
  }
  S <- svd(B)$d
  S <- sort(S)

  if (min(rkA, rkB) > 0) {
    for (i in 1:min(rkA, rkB)) {
      if (C[i]^2 < 0.5) {
        angles[1, i] <- acos(C[i])
      } else if (S[i]^2 <= 0.5) {
        angles[1, i] <- asin(S[i])
      }
    }
    angles <- t(angles)
  }
  list(angles = angles)
}


#' True positive rate
#'
#' Compares the sparsity structure of an estimated matrix `A` with the ground
#' truth `A_true`.
#'
#' @param A_true A matrix of the ground truth.
#' @param A A matrix of the same dimension.
#' @returns The proportion of non-zero entries of `A_true` that are also
#'  non-zero in `A`.
#' @export
#' @examples
#' A_true <- matrix(c(1, 1, 0, 0), 2)
#' A <- matrix(c(1, 0, 0, 1), 2)
#' TPR(A_true, A)
TPR <- function(A_true, A) {
  sum(A_true != 0 & A != 0) / sum(A_true != 0)
}


#' True negative rate
#'
#' Compares the sparsity structure of an estimated matrix `A` with the ground
#' truth `A_true`.
#'
#' @param A_true A matrix of the ground truth.
#' @param A A matrix of the same dimension.
#' @returns The proportion of zero entries of `A_true` that are also zero in
#'  `A`.
#' @export
#' @examples
#' A_true <- matrix(c(1, 1, 0, 0), 2)
#' A <- matrix(c(1, 0, 0, 1), 2)
#' TNR(A_true, A)
TNR <- function(A_true, A) {
  sum(A_true == 0 & A == 0) / sum(A_true == 0)
}

#' Higher order constraint for robust and sparse PCA / CCA
#'
#' Function to compute higher order constraint for a tensor
#'
#' @param z       A (1xp) torch::torch_tensor
#' @param low_z   A (kxp) torch::torch_tensor
#' @param C       A (pxp) matrix
#' @noRd
higher_order <- function(z, low_z, C = NULL) {

  if (is.null(C)) {
    low_u <- -torch::torch_sum(torch::torch_abs(torch::torch_matmul(torch::torch_transpose(z, 1, -1), low_z)), 2)
  } else{
    low_u <- -torch::torch_sum(torch::torch_abs(torch::torch_matmul(torch::torch_transpose(z, 1, -1), torch::torch_matmul(C, low_z))), 2)}
  return(torch::torch_unsqueeze(low_u, 2))
}

#' Domain constraint for robust and sparse PCA / CCA
#'
#' Function to compute domain constraint for a tensor
#'
#' @param z    A (1xp) torch::torch_tensor
#' @param C    A (pxp) quadratic matrix
#' @returns z_constrained, the value of the equality constraints
#' @noRd
constraint_fun <- function(z, C = NULL) {
  if (is.null(C)) {
    z_constrained <- torch::torch_matmul(torch::torch_transpose(z, 1, -1), z) - 1
  } else {
    z_constrained <- torch::torch_matmul(torch::torch_transpose(z, 1, -1), torch::torch_matmul(C, z)) - 1
  }
  if (as.numeric(z_constrained) > 0) {
    return(z_constrained)
  } else {
    return(torch::torch_unsqueeze(torch::torch_tensor(c(0)), 2))
  }
}

#' Penalty constraints for robust and sparse PCA / CCA
#'
#' Function to compute elastic net penalty for a tensor
#'
#' @param z             A torch::torch_tensor
#' @param pen           penalty parameter
#' @param alpha         elastic net mixing parameter
#' @returns              A 1x1 torch::torch_tensor
#' @noRd
penalty_fun <- function(z, pen, alpha) {
  l1 <- torch::torch_sum(torch::torch_abs(z))
  l2 <- torch::torch_square(torch::torch_norm(z))
  penalty <- alpha * l1 + (1 - alpha) * l2 - pen

  if (as.numeric(penalty) > 0) {
    return(torch::torch_unsqueeze(penalty, 2))
  } else {
    return(torch::torch_unsqueeze(torch::torch_tensor(c(0)), 2))
  }
}

#' Penalty constraints for robust and sparse PCA via SVD
#'
#' Function to compute elastic net penalty for a tensor
#'
#' @param z             A torch::torch_tensor
#' @param pen           penalty parameter
#' @param alpha         elastic net mixing parameter
#' @returns              A 1x1 torch::torch_tensor
#' @noRd
penalty_fun_unc <- function(z, pen, alpha) {
  l1 <- torch::torch_sum(torch::torch_tanh(1000*z) * z)
  l2 <- torch::torch_square(torch::torch_norm(z))
  penalty <- (alpha * l1 + (1 - alpha) * l2) * pen

  return(torch::torch_unsqueeze(penalty, 2))

}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute Pseudo-Huber loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' @param d             Scalar, parameter for Huber loss
#' @returns             A 1x1 torch::torch_tensor
#' @noRd
loss_huber <- function(M, d = 1.35) {
  sigma <- torch::torch_median(torch::torch_abs(M), 1, keepdim = TRUE)[[1]]
  M <- M / sigma
  loss <- torch::torch_mean(d^2 * (torch::torch_sqrt(1 + (M/d)^2) - 1) * sigma^2)
  return(loss)
}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute Tukey loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' @param k             Scalar, tuning constant of the Tukey loss
#' @returns             A 1x1 torch::torch_tensor
#' @noRd
loss_tukey <- function(M, k = 1.35) {
  sigma <- torch::torch_median(torch::torch_abs(M), 1, keepdim = TRUE)[[1]]
  M <- M / sigma
  tukey <- (M/k)^2 * (3 - 3*(M/k)^2 + (M/k)^4)
  loss <- torch::torch_mean(torch::torch_where(torch::torch_abs(M) <= k, tukey, 1.) * sigma^2)
  return(loss)
}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute trimmed loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' Per column, only the proportion `h` of cells with the smallest squared
#' residuals contributes to the loss; the other cells are trimmed.
#'
#' @param h             Trimming parameter (robustness) of observations to use per column
#' @returns             A 1x1 torch::torch_tensor
#' @noRd
loss_lts <- function(M, h = 0.75) {
  least_squares <- M^2
  # h-quantile of the squared residuals of each column; before version 0.2.0
  # the quantile of the signed residuals M was used and compared with `<`,
  # which kept far fewer cells than intended
  threshold <- torch::torch_quantile(least_squares,
                                     torch::torch_tensor(h, dtype = torch::torch_double()),
                                     dim = 1,
                                     keepdim = TRUE,
                                     interpolation = "nearest")
  lts_loss <- torch::torch_where(least_squares <= threshold, least_squares, 0)
  loss <- torch::torch_mean(lts_loss)
  return(loss)
}


#' F1 score
#'
#' Function to compute F1 score for outlier detection, given the true positions
#'
#' @param ind_outl     position of outliers
#' @param flag_outl    flagged cells
#' @noRd
f1_score <- function(ind_outl, flag_outl) {
  TP <- length(intersect(flag_outl, ind_outl)) #TP
  FP <- length(setdiff(flag_outl, ind_outl)) #FP
  FN <- length(setdiff(ind_outl, flag_outl)) #FN

  TP / (TP + 0.5 * (FP + FN))
}
