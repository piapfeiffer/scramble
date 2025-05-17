# Hyperparameter Optimization for routines in scramble

# ==============================================
# Functions for PCA
# ==============================================
#' Bayesian Hyperparameter optimization for PCA via `scramble` algorithm
#'
#' Get optimal penalties using bayesian optimization.
#' Score function depends on whether sparsity is assumed or not
#' @param Q Q from QR decomposition on data_x
#' @param R R from QR decomposition on data_x
#' @param rsvd object resulting from SVD on R
#' @param p number of variables in data_x
#' @param n number of observations
#' @param alpha_x scalar: A positive real scalar between 0 and
#' 1, indicating the elastic
#'  net parameter for data_x.
#' @param tol tolerance for training algorithm
#' @param lr  learning rate for gradient descent training algorithm
#' @param epochs maximum number of epochs for gradient descent training algorithm
#' @param lr_decay parameter for
#' @param criterion type of criterion used to determine optimal hyperparameters,
#' default: criterion = TPO
#' @returns a list object containing the optimal hyperparameters and a summary
#' of the optimization process
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
                              0,
                              0,
                              loss_type, h,
                              tol = tol, lr = lr, epochs = epochs,
                              lr_decay = lr_decay)}

  scoringFunction <- function(pen_x){
    # PEN_X <- c(...)
    res_train <- train_PCA_SVD(X, groups,
                               rsvd, n, p,
                               pen_x,
                               alpha_x,
                               loss_type, h,
                               tol = tol, lr = lr, epochs = epochs,
                               lr_decay = lr_decay)

    R <- res_train$best_measure
    tot_var <- sum(apply(X, 2, robustbase::Qn))
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
  bounds <- list()

  bounds <- list(pen_x = bounds_input)

  summary <- c()
  optObj <- c()
  optObj$scoreSummary <- NA

  best_params <- list()

  best_params$pen_x <- list(pen_x = max(bounds$pen_x))
  summary_list <- list()

  out <- tryCatch({

    optObj <- ParBayesianOptimization::bayesOpt(
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
    )
    best_params <- ParBayesianOptimization::getBestPars(optObj)
  },  error = function(err) {
    print(paste("ERROR:  ",err, "\n maximum penalty values are used!"))
  }, finally = {
    return(list(best_params = best_params, summary = optObj$scoreSummary, proc = optObj$GauProList))
  })
}
