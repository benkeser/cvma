#' Helper function to get validation fold predictions. Each
#' element of the \code{fit$pred} list will have a name
#' corresponding to the validation fold the predictions
#' were obtained on. 
#' @keywords internal 
get_valid_pred_from_fit <- function(fit, valid_folds){
    unlist(fit$pred[as.character(valid_folds)])
}

#' Helper function to get validation fold outcomes. 
#' @keywords internal
get_valid_y <- function(split_Y, valid_folds){
    split_Y[as.character(valid_folds)]
}

#' Helper function to get a vector of fold assignments (deprecated)
#' @keywords internal
get_fold_vector <- function(valid_y){
    if(is.null(dim(valid_y[[1]]))){
        as.numeric(rep(names(valid_y), lapply(valid_y, length)))        
    }else{
        as.numeric(rep(names(valid_y), lapply(valid_y, nrow)))
    }
}

#' Helper function to get outcomes in proper format needed
#' for some of the \code{get_} functions. 
#' @keywords internal 
get_Y_out <- function(x, split_Y, training_folds, V = NULL){
    if(!is.null(training_folds)){
        valid_folds <- training_folds[-which(training_folds %in% x)]        
    }else{
        valid_folds <- (1:V)[-x]
    }
    valid_y <- get_valid_y(split_Y, valid_folds)
    return(Reduce(rbind, valid_y))
}

#' Helper function to format the validation folds used in
#' some of the \code{get_} functions
#' @keywords internal 
get_fold_out <- function(x, split_Y, training_folds, valid_folds = NULL, V = NULL){
    if(!is.null(x)){
        if(!is.null(training_folds)){
            # evaluates in get_sl_input
            # here, training_folds corresponds to the folds used to train 
            # the super learner, and we are getting which one of those folds
            # was held out in fitting each of the candidate learners  
            valid_folds <- training_folds[-which(training_folds %in% x)]
        }else{
            valid_folds <- (1:V)[-x]
        }
    }
    return(valid_folds)
}

#' Helper function to get learner prediction matrices formatted
#' properly. 
#' @keywords internal 
get_pred_out <- function(fit_idx, learner_folds = NULL, sl_folds = NULL, 
                         all_fits, valid_folds = NULL, learners){
    if(!is.null(learner_folds)){
        # evaluates in get_sl_input
        # where we want to identify validation folds for the super learner 
        valid_folds <- sl_folds[-which(sl_folds %in% learner_folds)]            
    }
    Z <- Reduce(cbind, lapply(all_fits[fit_idx], get_valid_pred_from_fit, 
                              valid_folds = valid_folds))
    colnames(Z) <- learners
    return(Z)
}

#' Helper function to get super learner predictions formatted properly.
#' @keywords internal 
get_sl_pred_out <- function(sl_training_folds, Ynames, outer_valid_folds, 
                            all_fit_tasks, all_fits, all_sl, learners, 
                            sl_control, V){
  
    if(!is.null(outer_valid_folds)){
        # evaluates in get_y_weight_input
        # where we want to get sl predictions on the fold that was held out of 
        # the fitting of this sl, to use to get weights
        inner_valid_folds <- (1:V)[-c(sl_training_folds, outer_valid_folds)]        
    }else{
        # evaluates in get_risk and get_risk_sl
        inner_valid_folds <- (1:V)[-sl_training_folds]
    }
  
    # find fits corresponding with train_matrix[,x]
    fit_idx <- sapply(Ynames, search_fits_for_training_folds, training_folds = sl_training_folds, fits = all_fit_tasks,
                      simplify = FALSE)
    sl_idx <- sapply(Ynames, search_fits_for_training_folds, training_folds = sl_training_folds, fits = all_sl,
                      simplify = FALSE)
    
    # get prediction matrix
    pred_matrices <- lapply(fit_idx, get_pred_out, valid_folds = inner_valid_folds,
                            all_fits = all_fits, learners = learners)

    sl_weight <- lapply(sl_idx, FUN = function(x){
        all_sl[[x]]$sl_weight
    })

    # get super learner prediction
    sl_pred_list <- mapply(p = pred_matrices, s = sl_weight, FUN = function(p,s){
        do.call(sl_control$ensemble_fn, args = list(weight = s, pred = p))
    }, SIMPLIFY = FALSE)
        
    sl_pred_out <- Reduce(cbind, sl_pred_list)

    return(sl_pred_out)
}

