# Hyperparameter Optimization for routines in scramble

# ==============================================
# Functions for PCA
# ==============================================
#' Bayesian hyperparameter optimization for PCA via the `scramble` algorithm
#'
#' Finds the penalty parameter `pen_x` by Bayesian optimization of a score that
#' depends on `criterion`. If the optimization fails, a warning is given and
#' the upper bound of `bounds_input` is used.
#'
#' @param X,groups,rsvd,n,p,alpha_x,loss_type,h,tol,lr,epochs,lr_decay Passed
#'  on to [train_PCA_SVD()].
#' @param k Number of components; determines the number of evaluations.
#' @param bounds_input Numeric vector of length 2, the search range of `pen_x`.
#' @param criterion `"TPO"` (tradeoff product optimization), `"BIC"`, or any
#'  other value to maximize the explained variance only.
#' @returns A list with the optimal hyperparameters (`best_params`), a summary
#'  of the optimization (`summary`) and the fitted Gaussian process (`proc`).
#' @noRd
bayesian_optimization_PCA_SVD <- function(X, groups,
                                          rsvd, n, p, k,
                                      alpha_x,
                                      loss_type, h,
                                      tol, lr, epochs,
                                      lr_decay,
                                      bounds_input,
                                      criterion = "TPO"
){
  get_sparsity <- function(a){return(sum(a!=0))}

  if (criterion == "BIC"){
    res_full <- train_PCA_SVD(X, groups,
                              rsvd, n, p,
                              pen_x = 0,
                              alpha_x = 0,
                              loss_type = loss_type, h = h,
                              tol = tol, lr = lr, epochs = epochs,
                              lr_decay = lr_decay)
  }

  scoringFunction <- function(pen_x){
    res_train <- train_PCA_SVD(X, groups,
                               rsvd, n, p,
                               pen_x = pen_x,
                               alpha_x = alpha_x,
                               loss_type = loss_type, h = h,
                               tol = tol, lr = lr, epochs = epochs,
                               lr_decay = lr_decay)

    R <- res_train$best_measure
    p_rate <- apply(res_train$best_v, 2, get_sparsity)

    if (criterion == "TPO"){
      ev <- abs(R)
      sparsity <- 1 - p_rate / p + 1e-5
      score <- sum(ev * sparsity)
    } else if (criterion == "BIC"){
      A_full <- res_full$best_v
      A_const <- res_train$best_v

      var_full <- sum(apply(X - X %*% (A_full %*% t(A_full)), 2, robustbase::Qn))
      var_constrained <- sum(apply(X - X %*% (A_const %*% t(A_const)), 2, robustbase::Qn))

      score <- var_constrained / var_full + sum(p_rate) * log(n)/n

    } else {
      score <- abs(R)
    }

    pred <- sum(abs(R))

    return(list(Score = score, Pred = pred,
                Sparse_x = sum(1-p_rate/p)
    ))
  }
  bounds <- list(pen_x = bounds_input)

  opt_obj <- tryCatch(
    ParBayesianOptimization::bayesOpt(
      FUN = scoringFunction
      , bounds = bounds
      , acq = "ei"
      , initPoints = max(c(k + 1, 5))
      , iters.n = max(c(2 * k, 15))
      , eps = 0.1
      , parallel = FALSE
      , errorHandling = 3
      , acqThresh = 1
      , plotProgress = FALSE
      , verbose = 0
    ),
    error = function(err) {
      warning("Hyperparameter optimization failed: ", conditionMessage(err),
              "\nThe maximum penalty value is used.", call. = FALSE)
      NULL
    })

  if (is.null(opt_obj)) {
    return(list(best_params = list(pen_x = max(bounds$pen_x)),
                summary = NA,
                proc = NULL))
  }
  list(best_params = ParBayesianOptimization::getBestPars(opt_obj),
       summary = opt_obj$scoreSummary,
       proc = opt_obj$GauProList)
}
