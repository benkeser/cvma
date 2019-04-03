#' Changing options for super learner and outcome weighting schemes after a call
#' to \code{cvma}. 
#' 
#' This function allows you to use the learner fits from a previous call to \code{cvma}
#' to change options for how the super learner weights and outcome weights are computed. 
#' The majority of the computation time in \code{cvma} is spent fitting learners, while 
#' the (re-)computation of weights is relatively quick. Thus, this function allows one
#' to obtain new results without the need to re-fit a trove of learners. 
#' 
#' @param object An object of class \code{cvma}. 
#' @param Y A matrix or data.frame of outcomes. This is assumed to be the same matrix
#' or data.frame of outcomes used in the original call to \code{cvma} and the function
#' will not check whether this is true. It is the responsibility of the user to ensure that
#' this function is invoked with the correct data set. 
#' @param X A matrix or data.frame of predictors. This is assumed to be the same matrix
#' or data.frame of outcomes used in the original call to \code{cvma} and the function
#' will not check whether this is true. It is the responsibility of the user to ensure that
#' this function is invoked with the correct data set.
#' @param sl_control A list with named entries ensemble_fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param y_weight_control A list with named entries ensemble_fn, optim_risk_fn, weight_fn,
#' cv_risk_fn. Available functions can be viewed with \code{y_weight_control_options()}. See
#' \code{?y_weight_control_options} for more on how users may supply their own functions.  
#' @param return_control A list with named entries \code{outer_weight} (whether to return outcome
#' weights for outer-most fold of CV, default \code{TRUE}), \code{outer_sl} (whether to return the 
#' super learner fit for each outcome on all the data), \code{all_y} (whether to return cross-validated
#' performance metrics for all outcomes), \code{all_learner_assoc} (whether to return cross-validation
#' performance metrics for all learners), \code{all_learner_fits} (whether to return all learner fits, 
#' which, while memory intensive, can be helpful if association measures based on different outcome 
#' weighting schemes are desired). 
#' @param scale Standardize each outcome to be mean zero with standard deviation 1. This is assumed 
#' to be the same value as in the original call to \code{cvma} and the function
#' will not check whether this is true. It is the responsibility of the user to ensure that
#' this function is invoked with the correct option to match the original call. 
#' 
#' @return If \code{return_outer_sl} is TRUE, it will return for each outcome Super Learner fit weights 
#' and associated risk for each learner. In addition, it will return the fit for all learners based on 
#' all folds. If \code{return_outer_weight} is TRUE, it will return the weights for each outcome
#' obtained using V-1 cross-validation. If \code{return_all_y} is TRUE, it will return for each 
#' outcome cross-validated measure (nonparametric R-squared or AUC), confidence interval and associated
#' p-value. 
#' 
#' TO DO: Re-write this documentation. 
#'  
#' @export
#' 
#' @seealso predict method 
#' 
#' @importFrom future.apply future_lapply
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
#' # results for super learner and R^2 for convex 
#' # combination of outcomes
#' fit <- cvma(Y = Y, X = X, V = 5, 
#'             learners = c("SL.glm","SL.mean"), 
#' 			   return_control = list(outer_weight = TRUE,
#'                                   outer_sl = TRUE,
#'                                   inner_sl = FALSE, 
#'                                   all_y = TRUE,
#'                                   all_learner_assoc = TRUE,
#'                                   all_learner_fits = TRUE))
#' 
#' # re-weight with discrete super learner and R^2 
#' # for a single outcome
#' re_fit <- reweight_cvma(fit, Y = Y, X = X, 
#' 				sl_control = list(ensemble_fn = "ensemble_linear",
#'                                optim_risk_fn = "optim_risk_sl_se",
#'                                weight_fn = "weight_sl_01",
#'                                cv_risk_fn = "cv_risk_sl_r2",
#'                                family = gaussian(),
#'                                alpha = 0.05),
#' 				y_weight_control = list(ensemble_fn = "ensemble_linear",
#'                                  weight_fn = "weight_y_01",
#'                                  optim_risk_fn = "optim_risk_y_r2",
#'                                  cv_risk_fn = "cv_risk_y_r2",
#'                                  alpha = 0.05))


reweight_cvma <- function(object, Y, X, 
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
                      return_control = list(outer_weight = TRUE,
                                            outer_sl = TRUE,
                                            inner_sl = FALSE, 
                                            all_y = TRUE,
                                            all_learner_assoc = TRUE,
                                            all_learner_fits = FALSE),
                      scale = FALSE
                      ){
    
    # get initial parameter values
    n <- length(Y[,1])
    J <- ncol(Y)
    M <- length(object$learners)
    
    #scale all multivariate outcome by covariance matrix?
    if(scale & !(all(Y==0 | Y==1))){
      Ymat <- data.frame(scale(Y))
    }else{
      # put Y into the proper format
      Ymat <- data.matrix(Y)
    }
    
    # Checks
    if(!all(apply(Ymat,2,function(x) { all(x %in% 0:1) })) & sl_control$optim_risk_fn == "optim_risk_sl_auc"){
      stop("Outcome should be binary.")
    }
    
    if(y_weight_control$weight_fn == "weight_y_convex" & y_weight_control$optim_risk_fn == "optim_risk_y_auc"){
      stop("risk_y_auc requires all composite outcome to be either 0 or 1.")
    }

    if(is.null(object$all_learner_fits)){
    	stop("reweight_cvma requires non-NULL object$all_learner_fits. Re-run cvma with return_control$return_learner_fits = TRUE.")
    }
    # correct names if none
    if(is.null(colnames(Ymat))){
        colnames(Ymat) <- paste0("Y",1:J)
    }
    
    folds <- object$folds
    all_fit_tasks <- make_fit_task_list(Ynames = colnames(Ymat), learners = object$learners, 
                                    V = object$V, return_outer_sl = return_control$outer_sl)
    
    # all super learner weight-getting tasks
    all_sl_tasks <- make_sl_task_list(Ynames = colnames(Ymat), V = object$V)
    # TO DO: I have a hunch that if future_lapply requires transferring
    #        all_fits between nodes that the communication overhead will make
    #        parallelization of this step slower than doing it sequentially 
    all_sl <- lapply(all_sl_tasks, FUN = get_sl, 
                            Y = Ymat, V = object$V, all_fit_tasks = all_fit_tasks, 
                            all_fits = object$all_learner_fits, folds = object$folds,
                            learners = object$learners, sl_control = sl_control)

    # all outcome weight tasks
    if(length(colnames(Ymat)) > 1){
      all_y_weight_tasks <- make_y_weight_task_list(V = object$V)      
    }else{
      all_y_weight_tasks <- list(list(training_folds = 1:object$V))
    }
    all_weight <- lapply(all_y_weight_tasks, FUN = get_y_weight, 
                                 Y = Ymat, V = object$V, Ynames = colnames(Ymat), 
                                 all_fits = object$all_learner_fits, all_sl = all_sl, 
                                 all_fit_tasks = all_fit_tasks, 
                                 sl_control = sl_control, y_weight_control = y_weight_control, 
                                 folds = object$folds,
                                 learners = object$learners)

    # get risk of entire procedure
    risk <- get_risk(Y = Y, V = object$V, all_fit_tasks = all_fit_tasks, 
                     all_fits = object$all_learner_fits, all_sl = all_sl, folds = object$folds, 
                     all_weight = all_weight, 
                     sl_control = sl_control, y_weight_control = y_weight_control, 
                     learners = object$learners) 

    # compute outer super learner
    if(return_control$outer_sl){
        outer_sl_tasks <- make_outer_sl_task_list(Ynames = colnames(Ymat), V = object$V)
        all_outer_sl <- lapply(outer_sl_tasks, FUN = get_formatted_sl, 
                           Y = Ymat, V = object$V, all_fit_tasks = all_fit_tasks,
                           all_fits = object$all_learner_fits, folds = object$folds,
                           sl_control = sl_control, learners = object$learners, 
                           return_learner_fits = TRUE)
    }else{
        all_outer_sl <- NULL
    }

    # compute outer weights
    if(return_control$outer_weight){
        outer_weight <- get_y_weight(task = list(training_folds = 1:object$V),
                                       Y = Ymat, V = object$V, Ynames = colnames(Ymat), 
                                       all_fits = object$all_learner_fits, all_sl = all_sl, 
                                       all_fit_tasks = all_fit_tasks,
                                       sl_control = sl_control, y_weight_control = y_weight_control, folds = object$folds,
                                       learners = object$learners)
    }else{
        outer_weight <- NULL
    }

    # get CV risk for each learner and each outcome 
    if(return_control$all_y){
      outer_sl_tasks <- make_outer_sl_task_list(Ynames = colnames(Ymat), V = object$V)
      risk_all_y <- lapply(outer_sl_tasks, FUN = get_risk_sl, 
                           Y = Y, V = object$V, all_fit_tasks = all_fit_tasks, 
                           all_fits = object$all_learner_fits, all_sl = all_sl,
                           folds = object$folds, sl_control = sl_control, learners = object$learners)
    }else{
        risk_all_y <- NULL
    }

    if(return_control$all_learner_assoc){
      outer_learner_tasks <- make_outer_learner_task_list(Ynames = colnames(Ymat),
                                                          V = object$V, learners = object$learners)
      risk_all_learners <- lapply(outer_learner_tasks, FUN = get_risk_learner, 
                           Y = Y, V = object$V, all_fit_tasks = all_fit_tasks, 
                           all_fits = object$all_learner_fits, all_sl = all_sl,
                           folds = object$folds, sl_control = sl_control, learners = object$learners)
    }else{
      risk_all_learners <- NULL
    }

    if(!return_control$all_learner_fits){
      all_fits <- NULL
    }

    # format output
    out <- list(cv_assoc = risk, 
                sl_fits = all_outer_sl, 
                outer_weight = outer_weight,
                inner_weight = all_weight, 
                cv_assoc_all_y = risk_all_y,
                cv_assoc_all_learners = risk_all_learners,
                folds = object$folds, 
                y_names = colnames(Ymat),
                learners = object$learners,
                all_learner_fits = all_fits)
    class(out) <- "cvma"
    return(out)
} 