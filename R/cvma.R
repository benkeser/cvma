#' Cross-validated maximal association measures
#' 
#' A flexible interface for computing cross-validation-based measures of maximal association.
#' In an outer layer of V-fold cross validation, training samples are used to train a prediction 
#' algorithm for each outcome. Multiple algorithms may be ensembled using stacking (also known as 
#' super learning) based on V-2 fold cross-validation. An inner layer of V-1 cross validation is used 
#' to determine a user-specified combination of outcomes that maximizes a user-specified prediction 
#' criteria. The outer layer
#' validation sample is used to compute a user-specified cross-validated measure of performance of
#' the prediction algorithm for predicting the combined outcome that was computed in the training 
#' sample. Several common choices for outcome combinations (convex combination of
#' outcomes and single outcome that is most associated) and prediction criteria (nonparametric R^2,
#' negative log-likelihood, and area under ROC curve) are included; however, users may specify
#' their own criteria as well. The function returns the cross-validated summary measure for the
#' maximally combined outcome and, if desired, the cross-validated summary measure for each 
#' outcome.  
#' 
#' TO DO: Figure out how future works (e.g., can plan() be specified internally
#' or externally?)
#' 
#' @param Y A matrix or data.frame of outcomes
#' @param X A matrix or data.frame of predictors
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param learners Super learner wrappers. See \code{SuperLearner::listWrappers}.
#' @param sl_control A list with named entries ensemble_fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param y_weight_control A list with named entries ensemble_fn, optim_risk_fn, weight_fn,
#' cv_risk_fn. Available functions can be viewed with \code{y_weight_control_options()}. See
#' \code{?y_weight_control_options} for more on how users may supply their own functions.  
#' @param return_outer_weight Whether to return estimate of outcome
#' weights the outer most cross-validation layer (i.e, based on V - 1 cross-validated 
#' super learner risks).
#' @param return_outer_sl Whether to return the super learner for the outer most cross-validation
#' layer (i.e., V-fold super learner).
#' @param return_all_y Whether to return cross-validated performance measures for each
#' column of \code{Y}. 
#' @param scale Standardize each outcome to be mean zero with standard deviation 1.
#' 
#' @return If \code{return_outer_sl} is TRUE, it will return for each outcome Super Learner fit weights 
#' and associated risk for each learner. In addition, it will return the fit for all learners based on 
#' all folds. If \code{return_outer_weight} is TRUE, it will return the weights for each outcome
#' obtained using V-1 cross-validation. If \code{return_all_y} is TRUE, it will return for each 
#' outcome cross-validated measure (nonparametric R-squared or AUC), confidence interval and associated
#' p-value. 
#'  
#' @export
#' 
#' @seealso predict method 
#' 
#' @importFrom future future_lapply
#' @importFrom stats gaussian
#' @import SuperLearner 
#' 
#' @examples
#' set.seed(1234)
#' library(SuperLearner)
#' library(future)
#' X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
#' Y1 <- rnorm(100, X$x1 + X$x2, 1)
#' Y2 <- rnorm(100, X$x1 + X$x2, 3)
#' Y <- data.frame(Y1 = Y1, Y2 = Y2)
#' fit <- cvma(Y = Y, X = X, V = 5, 
#'                 learners = c("SL.glm","SL.mean"))
#' 
#' 
cvma <- function(Y, X, V = 5, learners, 
                      sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_se",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_r2",
                                   family = gaussian(),
                                   alpha = 0.05),
                      y_weight_control = list(ensemble_fn = "ensemble_linear",
                                  weight_fn = "weight_y_convex",
                                  optim_risk_fn = "optim_risk_y_r2",
                                  cv_risk_fn = "cv_risk_y_r2",
                                  alpha = 0.05),
                      return_outer_weight = TRUE,
                      return_outer_sl = TRUE,
                      return_all_y = TRUE,
                      scale=FALSE
                      ){
    
    # get initial parameter values
    n <- length(Y[,1])
    J <- ncol(Y)
    M <- length(learners)
    
    #scale all multivariate outcome by covariance matrix?
    if(scale & !(all(Y==0 | Y==1))){
      Ymat <- data.frame(scale(Y))
    }else{
      # put Y into the proper format
      Ymat <- data.matrix(Y)
    }
    
    #Checks
    if(!all(apply(Ymat,2,function(x) { all(x %in% 0:1) })) & sl_control$optim_risk_fn == "optim_risk_sl_auc"){
      stop("Outcome should be binary.")
    }
    
    if(y_weight_control$weight_fn == "weight_y_convex" & y_weight_control$optim_risk_fn == "optim_risk_y_auc"){
      stop("Risk_y_auc requires all composite outcome to be either 0 or 1.")
    }

    # correct names if none
    if(is.null(colnames(Ymat))){
        colnames(Ymat) <- paste0("Y",1:J)
    }
    
    # TO DO: make_folds function with more options, possibly from origami?
    folds <- rep(seq_len(V), length = n)
    folds <- sample(folds)

    # all learner fitting tasks
    all_fit_tasks <- make_fit_task_list(Ynames = colnames(Ymat), learners = learners, 
                                        V = V, return_outer_sl = return_outer_sl)
    
    all_fits <- future::future_lapply(all_fit_tasks, FUN = get_fit, folds = folds, 
                              X = X, Y = Ymat, sl_control = sl_control)

    # all super learner weight-getting tasks
    all_sl_tasks <- make_sl_task_list(Ynames = colnames(Ymat), V = V)
    # TO DO: I have a hunch that if future_lapply requires transferring
    #        all_fits between nodes that the communication overhead will make
    #        parallelization of this step slower than doing it sequentially 
    all_sl <- lapply(all_sl_tasks, FUN = get_sl, 
                            Y = Ymat, V = V, all_fit_tasks = all_fit_tasks, 
                            all_fits = all_fits, folds = folds,
                            learners = learners, sl_control = sl_control)

    # all outcome weight tasks
    all_y_weight_tasks <- make_y_weight_task_list(V = V)
    all_weight <- lapply(all_y_weight_tasks, FUN = get_y_weight, 
                                 Y = Ymat, V = V, Ynames = colnames(Ymat), 
                                 all_fits = all_fits, all_sl = all_sl, 
                                 all_fit_tasks = all_fit_tasks, 
                                 sl_control = sl_control, y_weight_control = y_weight_control, folds = folds,
                                 learners = learners)

    # get risk input
    risk <- get_risk(Y = Y, V = V, all_fit_tasks = all_fit_tasks, 
                     all_fits = all_fits, all_sl = all_sl, folds = folds, 
                     all_weight = all_weight, 
                     sl_control = sl_control, y_weight_control = y_weight_control, learners = learners) 

    # compute outer super learner
    if(return_outer_sl){
        outer_sl_tasks <- make_outer_sl_task_list(Ynames = colnames(Ymat), V = V)
        all_outer_sl <- lapply(outer_sl_tasks, FUN = get_formatted_sl, 
                           Y = Ymat, V = V, all_fit_tasks = all_fit_tasks,
                           all_fits = all_fits, folds = folds,
                           sl_control = sl_control, learners = learners, 
                           return_learner_fits = TRUE)
    }else{
        all_outer_sl <- NULL
    }

    # compute outer weights
    if(return_outer_weight){
        outer_weight <- get_y_weight(task = list(training_folds = 1:V),
                                       Y = Ymat, V = V, Ynames = colnames(Ymat), 
                                       all_fits = all_fits, all_sl = all_sl, 
                                       all_fit_tasks = all_fit_tasks,
                                       sl_control = sl_control, y_weight_control = y_weight_control, folds = folds,
                                       learners = learners)
    }else{
        outer_weight <- NULL
    }

    # get CV risk for each outcome 
    if(return_all_y){
        outer_sl_tasks <- make_outer_sl_task_list(Ynames = colnames(Ymat), V = V)
        risk_all_y <- lapply(outer_sl_tasks, FUN = get_risk_sl, 
                             Y = Y, V = V, all_fit_tasks = all_fit_tasks, 
                             all_fits = all_fits, all_sl = all_sl,
                             folds = folds, sl_control = sl_control, learners = learners)
    }else{
        risk_all_y <- NULL
    }

    # format output
    out <- list(cv_assoc = risk, 
                sl_fits = all_outer_sl, 
                outer_weight = outer_weight,
                inner_weight = all_weight, 
                cv_assoc_all_y = risk_all_y,
                y_names = colnames(Ymat))
    class(out) <- "cvma"
    return(out)
} 