#' Cross-validated maximal association measures
#' 
#' TO DO: Add a description here.
#' 
#'  
#'    
#' TO DO: Figure out how future works (e.g., can plan() be specified internally
#' or externally?)
#' 
#' @param Y A matrix or data.frame of outcomes
#' @param X A matrix or data.frame of predictors
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param learners Super learner wrappers. See \code{SuperLearner::listWrappers}.
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param y_weight_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
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
#' @return TO DO : Add return values
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
    if(scale){
      Ymat <- data.frame(scale(Y))
    }else{
      # put Y into the proper format
      Ymat <- data.matrix(Y)
    }
    
    #Checks
    if(!all(apply(Ymat,2,function(x) { all(x %in% 0:1) })) & sl_control$optim_risk_fn == "optim_risk_sl_auc"){
      stop("Outcome should be binary.")
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
                                 # ensemble_fn = ensemble_fn,
                                 # optim_risk_y_weight_control = optim_risk_y_weight_control)

    # get risk input
    risk <- get_risk(Y = Y, V = V, all_fit_tasks = all_fit_tasks, 
                     all_fits = all_fits, all_sl = all_sl, folds = folds, 
                     all_weight = all_weight, 
                     sl_control = sl_control, y_weight_control = y_weight_control, learners = learners) 
                     # ensemble_fn = ensemble_fn, cv_risk_y_weight_control = cv_risk_y_weight_control)

    # compute outer super learner
    if(return_outer_sl){
        outer_sl_tasks <- make_outer_sl_task_list(Ynames = colnames(Ymat), V = V)
        all_outer_sl <- lapply(outer_sl_tasks, FUN = get_formatted_sl, 
                           Y = Ymat, V = V, all_fit_tasks = all_fit_tasks,
                           all_fits = all_fits, folds = folds,
                           sl_control = sl_control, learners = learners, 
                           # ensemble_fn = ensemble_fn, risk_sl_control = risk_sl_control,
                           # weight_sl_control = weight_sl_control, 
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
                                       # ensemble_fn = ensemble_fn,
                                       # optim_risk_y_weight_control = optim_risk_y_weight_control)
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
                cv_assoc_all_y = risk_all_y)
    class(out) <- "max_assoc"
    return(out)
} 