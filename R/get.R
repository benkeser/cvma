#' Fit a learner on training folds and get predictions on validation folds
#' 
#' @param task A named list identifying what training_folds to fit the learner
#' on. The function returns predictions from this fit on the remaining folds (i.e.,
#' the validation folds). 
#' @param folds Vector identifying which fold observations fall into. 
#' @param Y A matrix or data.frame of outcomes.
#' @param X A matrix or data.frame of predictors.
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' 
#' @return A named list with task and output of super learner wrapper fit. 

get_fit <- function(task, folds, X, Y, sl_control){
    train_idx <- folds %in% task$training_folds
    valid_idx <- !train_idx

    train_Y <- Y[train_idx, task$Yname]
    train_X <- X[train_idx, , drop = FALSE]
    if(sum(valid_idx) > 0){
      valid_X <- X[valid_idx, , drop = FALSE]
    }else{
      # if learner being fit on all data, then return
      # predictions on training sample
      valid_X <- train_X
    }
    # fit super learner wrapper
    fit <- do.call(task$SL_wrap, list(Y = train_Y, X = train_X,
                                 newX = valid_X, obsWeights = rep(1, length(train_Y)),
                                 family = sl_control$family))
    # split up validation predictions
    if(sum(valid_idx) > 0){
      fit$pred <- split(fit$pred, folds[valid_idx])
    }else{
      fit$pred <- fit$pred
    }
    # return(fit)
    return(c(task, fit))
}

#' Get outcome weights based on cross-validated super learner fits
#' 
#' @param task A named list identifying what training folds to use to 
#' obtain outcome weights. 
#' @param Y A matrix or data.frame of outcomes
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param Ynames Names of the columns of \code{Y}.
#' @param all_fits A list of all learner fits (from \code{get_fit})
#' @param all_sl A list of all super learner fits (from \code{get_sl})
#' @param all_fit_tasks A list of all learner fitting tasks (quicker to search over
#' than \code{all_fits}). 
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param y_weight_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn. Available functions can be viewed with \code{y_weight_control_options()}. See
#' \code{?y_weight_control_options} for more on how users may supply their own functions.  
#' @param folds Vector identifying which fold observations fall into. 
#' @param learners Super learner wrappers. See \code{\link{SuperLearner::listWrappers}}.
#' 
#' @return Named list identifying training folds used and the composite outcome weights. 

get_y_weight <- function(task, Y, V, Ynames, all_fits, all_sl, 
                          all_fit_tasks, sl_control, y_weight_control, folds, learners){
    # add folds so split.data.frame works
    split_Y <- split(data.frame(Y, folds), folds) 
    # remove folds column
    split_Y <- lapply(split_Y, function(x){ x[,-ncol(x)] })

    if(!all(1:V %in% task$training_folds)){
        valid_folds <- (1:V)[-task$training_folds]
    }else{
        # evaluates when computing outer weights
        valid_folds <- NULL
    }
    # get input needed to compute weights
    input <- get_y_weight_input(valid_folds = valid_folds, 
                                 split_Y = split_Y, Ynames = Ynames, 
                                 all_fits = all_fits, all_sl = all_sl, 
                                 all_fit_tasks = all_fit_tasks, V = V, 
                                 sl_control = sl_control,
                                 learners = learners)

    y_weight <- do.call(y_weight_control$weight_fn, 
                         args = list(input = input, y_weight_control = y_weight_control))
    out <- c(list(training_folds = task$training_folds), y_weight)
    return(out)
}

#' Get super learner weights based on cross-validated learner fits
#' 
#' @param task A named list identifying what training folds to use to 
#' obtain outcome weights. 
#' @param Y A matrix or data.frame of outcomes
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param all_fits A list of all learner fits (from \code{get_fit})
#' @param all_fit_tasks A list of all learner fitting tasks (quicker to search over
#' than \code{all_fits}). 
#' @param sl_control A list with named entries ensemble_fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param folds Vector identifying which fold observations fall into. 
#' @param learners Super learner wrappers. See \code{\link{SuperLearner::listWrappers}}.
#' 
#' @return Named list identifying training folds used and the composite outcome weights. 
#' 
get_sl <- function(task, Y, V, all_fit_tasks, all_fits, folds, sl_control, learners){
    split_Y <- split(Y[ , task$Yname], folds)
    if(!all(1:V %in% task$training_folds)){
        valid_folds <- (1:V)[-task$training_folds]        
    }else{
        valid_folds <- NULL
    }
    # get input needed to compute sl ensemble weights 
    input <- get_sl_input(split_Y = split_Y, valid_folds = valid_folds,
                          Yname = task$Yname, V = V, all_fit_tasks = all_fit_tasks, 
                          all_fits = all_fits, learners = learners)
    # get sl ensemble weights
    sl_weight <- do.call(sl_control$weight_fn, 
                          args = list(input = input, sl_control = sl_control))

    out <- list(training_folds = task$training_folds, Yname = task$Yname, sl_weight = sl_weight$weight)
    # get sl predictions on valid_folds folds by searching for when
    # (1:V)[-valid_folds] is training sample
    # sl_pred <- get_sl_pred(valid_folds = valid_folds, V = V, all_fit_tasks = all_fit_tasks,
    #                        all_fits = all_fits, y = Y_name, sl_weight = sl_weight,
    #                        split_Y = split_Y)
    return(out)
}

#' Get a super learner fit for a given outcome with more information. 
#' 
#' This is called if \code{return_outer_sl = TRUE}, in which case more 
#' information on learner risks etc... is computed and returned. 
#' 
#' @param task A named list identifying what outcome to use. 
#' @param Y A matrix or data.frame of outcomes
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param all_fits A list of all learner fits (from \code{get_fit})
#' @param all_sl A list of all super learner fits (from \code{get_sl})
#' @param all_fit_tasks A list of all learner fitting tasks (quicker to search over
#' than \code{all_fits}). 
#' @param return_learner_fits Should the \code{fit} component of the learners
#' be returned. Must be \code{TRUE} to obtain later predictions. 
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param folds Vector identifying which fold observations fall into. 
#' @param learners Super learner wrappers. See \code{\link{SuperLearner::listWrappers}}.
#' @return Named list of super learner results. 

get_formatted_sl <- function(task, Y, V, all_fit_tasks, all_fits, folds, 
                             sl_control, return_learner_fits = TRUE, 
                             learners){
    split_Y <- split(Y[ , task$Yname], folds)
    if(!all(1:V %in% task$training_folds)){
        valid_folds <- (1:V)[-task$training_folds]        
    }else{
        valid_folds <- NULL
    }
    # get input needed to compute sl ensemble weights 
    input <- get_sl_input(split_Y = split_Y, valid_folds = valid_folds,
                          Yname = task$Yname, V = V, all_fit_tasks = all_fit_tasks, 
                          all_fits = all_fits, learners = learners)
    # get sl ensemble weights
    sl_weight_list <- do.call(sl_control$weight_fn, 
                          args = list(input = input, sl_control = sl_control))
    sl_weight <- sl_weight_list$weight

    # get risks for each learner
    all_pred <- Reduce(rbind, lapply(input, '[[', "pred"))
    all_y <- Reduce(c, lapply(input, '[[', "Y"))
    risk_input <- list(Y = all_y, pred = all_pred)

    M <- dim(input[[1]]$pred)[2]
    risks <- rep(Inf, M)
    for(m in 1:M){
        weights <- rep(0, M)
        weights[m] <- 1
        risks[m] <- do.call(sl_control$optim_risk_fn, 
                            args = list(sl_weight = weights, sl_control = sl_control,
                                        input = risk_input))
    }

    # get fits for each learner
    idx <- search_fits_for_training_folds(fits = all_fit_tasks, 
                                          y = task$Yname, 
                                          training_folds = 1:V)
    fit_out <- sapply(idx, function(i){
                lapply(all_fits[i], "[", "fit")
            })

    out <- list(Yname = task$Yname, sl_weight = sl_weight,
                learner_risks = risks, learers = fit_out)
    return(out)
}

#' Get cross-validated risk of the super learner for a univariate outcome
#' 
#' @param task A named list identifying what training folds to use to 
#' obtain outcome weights. 
#' @param Y A matrix or data.frame of outcomes
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param all_fits A list of all learner fits (from \code{get_fit})
#' @param all_sl A list of all super learner fits (from \code{get_sl})
#' @param all_fit_tasks A list of all learner fitting tasks (quicker to search over
#' than \code{all_fits}). 
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param folds Vector identifying which fold observations fall into. 
#' @param learners Super learner wrappers. See \code{\link{SuperLearner::listWrappers}}.
#' @return Named list identifying which outcome and the cross-validated risk of the super learner. 

get_risk_sl <- function(task, Y, V, all_sl, all_fit_tasks, all_fits, folds, 
                        sl_control, learners){
    split_Y <- split(Y[ , task$Yname], folds)
    
    if(!all(1:V %in% task$training_folds)){
        valid_folds <- (1:V)[-task$training_folds]        
    }else{
        valid_folds <- NULL
    }    
    input <- get_risk_sl_input(split_Y = split_Y, Yname = task$Yname, V = V,
                               all_fits = all_fits, 
                               all_weight = all_weight, all_sl = all_sl, 
                               all_fit_tasks = all_fit_tasks, folds = folds, 
                               sl_control = sl_control, learners = learners)
    risk <- do.call(sl_control$cv_risk_fn, args = list(input = input, sl_control = sl_control))
    out <- c(Yname = task$Yname, risk)
    return(out)
}

#' Get cross-validated risk of entire procedure (i.e., outer-most cross-validation layer)
#'  
#' @param Y A matrix or data.frame of outcomes
#' @param V Number of outer folds of cross-validation (nested cross-validation
#' uses V-1 and V-2 folds), so must be at least four. 
#' @param all_fits A list of all learner fits (from \code{get_fit})
#' @param all_sl A list of all super learner fits (from \code{get_sl})
#' @param all_fit_tasks A list of all learner fitting tasks (quicker to search over
#' than \code{all_fits}). 
#' @param all_weight A list of all outcome weights (from \code{get_y_weight})
#' @param sl_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn, family. Available functions can be viewed with \code{sl_control_options()}. See
#' \code{?sl_control_options} for more on how users may supply their own functions.  
#' @param y_weight_control A list with named entries ensemble.fn, optim_risk_fn, weight_fn,
#' cv_risk_fn. Available functions can be viewed with \code{y_weight_control_options()}. See
#' \code{?y_weight_control_options} for more on how users may supply their own functions.  
#' @param folds Vector identifying which fold observations fall into. 
#' @param learners Super learner wrappers. See \code{\link{SuperLearner::listWrappers}}.
#' 
#' @return Numeric risk value 

get_risk <- function(Y, V, all_fit_tasks, all_fits, all_weight, folds, all_sl,
                     sl_control, y_weight_control, learners){
    split_Y <- split(data.frame(Y, folds), folds) 
    split_Y <- lapply(split_Y, function(x){ x[,-ncol(x)] })
    input <- get_risk_input(split_Y = split_Y, Ynames = colnames(Y), 
                            all_fits = all_fits, all_weight = all_weight,
                            all_sl = all_sl, all_fit_tasks = all_fit_tasks, 
                            sl_control = sl_control, V = V,
                            learners = learners)
    risk <- do.call(y_weight_control$cv_risk_fn, 
                    args = list(input = input, y_weight_control = y_weight_control))
    return(risk)
}