#' Create input list for \code{get_y_weight}
#' 
#' @param split_Y The outcome matrix split by relevant validation folds.
#' @param Ynames The names of the outcomes. Used to search \code{all_fits}
#' and \code{all_sl}.
#' @param valid_folds If not \code{NULL}, identifies validation folds used in the
#' inner most super learner. Set to \code{NULL} when getting weights for outer
#' cross-validation folds (where all folds are used for training super learner).
#' @param all_fits List of all learner fits.
#' @param all_sl List of all super learner weight fits
#' @param all_fit_tasks List of all learner fit tasks (faster to search over than
#' search over all_fits). 
#' @param sl_control List of super learner controls.
#' @param V Number of folds.
#' @param learners Vector of super learner vectors.
#' 
#' @return List with each entry a list with entries: valid_fold (the number of the
#' corresponding fold), Y (the outcome matrix in this validation fold), pred (matrix
#' of outcome predictions for this validation fold from the super learner fit in 
#' training folds). The length of the list varies according to what is input into 
#' \code{valid_folds} so that this function may be used both in the outer and inner
#' layers of cross-validation. 

get_y_weight_input <- function(split_Y, Ynames, 
                                valid_folds, 
                                all_fits, all_sl, all_fit_tasks, 
                                sl_control, V, learners){
    # gives the splits that are used in the V-2 super learner
    if(!is.null(valid_folds)){
        training_folds <- (1:V)[-valid_folds]         
    }else{
        # evaluates when getting outer weights
        training_folds <- 1:V
    }
    # then we need to get fits for columns of this guy
    train_matrix <- utils::combn(training_folds, V - (length(valid_folds) + 1))
    # valid_folds <- training_folds[-which(training_folds %in% training_folds)]
    # need to search_fits_for_training_folds, but now for all Y
    # also need to search_sl_for_training_folds
    
    # Y output 
    Y_out <- lapply(split(train_matrix, col(train_matrix)), get_Y_out, 
                         split_Y = split_Y, training_folds = training_folds)

    # fold output 
    fold_out <- lapply(split(train_matrix, col(train_matrix)), 
                       get_fold_out, split_Y = split_Y, training_folds = training_folds,
                       valid_folds = valid_folds)

    # sl prediction for each outcome
    pred_out <- lapply(split(train_matrix, col(train_matrix)), 
                       get_sl_pred_out, V = V, 
                       Ynames = Ynames, outer_valid_folds = valid_folds, 
                       all_fit_tasks = all_fit_tasks, all_fits = all_fits, 
                       all_sl = all_sl, sl_control = sl_control, 
                       learners = learners)

    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_out, Y_out, pred_out)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred"); o})

    return(out)
}

#' Create input list for \code{get_sl}
#' 
#' @param split_Y The outcome matrix split by relevant validation folds.
#' @param Yname The names of the outcomes. Used to search \code{all_fits}
#' and \code{all_sl}.
#' @param valid_folds If not \code{NULL}, identifies validation folds used in the
#' inner most super learner. Set to \code{NULL} when getting weights for outer
#' cross-validation folds (where all folds are used for training super learner).
#' @param all_fits List of all learner fits.
#' @param all_fit_tasks List of all learner fit tasks (faster to search over than
#' search over all_fits). 
#' @param V Number of folds.
#' @param learners Vector of super learner vectors.
#' 
#' @return List with each entry a list with entries: valid_fold (the number of the
#' corresponding fold), Y (univariate outcome in this validation fold), pred (matrix
#' of outcome predictions for this validation fold from learners fit in 
#' training folds). The length of the list varies according to what is input into 
#' \code{valid_folds} so that this function may be used both in the outer and inner
#' layers of cross-validation. 

# this will loop over valid_folds values in columns of combn(V, 2)
get_sl_input <- function(split_Y, Yname, valid_folds, 
                         V, all_fit_tasks, all_fits, learners){
    # gives all splits used to compute this super learner
    if(!is.null(valid_folds)){
        sl_training_folds <- (1:V)[-valid_folds] 
    }else{
        # evaluates when getting outer super learner
        sl_training_folds <- 1:V
    }
    # gives a matrix whose columns correspond to the folds used 
    # to train each of the candidate 
    learner_training_folds <- utils::combn(sl_training_folds, V - (length(valid_folds) + 1))

    # find indexes corresponding to these training sets
    # each column corresponds to the indexes where the corresponding
    # column of train_matrix were in the training folds
    # each row corresponds to a different candidate learner 
    # TO DO: better book keeping could allow this to be done without
    #        using search_fits
    fit_idx_list <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                             search_fits_for_training_folds,
                             fits = all_fit_tasks, y = Yname)

    # let's assume that the method for obtaining weights will want an input
    # of Y, Y_hat (matrix from learners), and fold

    # need to put together matrix of predictions from learners.
    
    #Get actual outcome sorted by validation folds.
    Y_out <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                    get_Y_out, split_Y = split_Y, 
                    training_folds = sl_training_folds)

    # here x does get input to get_fold_out
    # valid_folds = folds in training_folds not used by the current fit
    fold_out <- lapply(split(learner_training_folds, col(learner_training_folds)), 
                       get_fold_out, split_Y = split_Y, 
                       training_folds = sl_training_folds)

    pred_out <- mapply(learner_folds = split(learner_training_folds, col(learner_training_folds)), 
                       fit_idx = fit_idx_list,
                       get_pred_out, 
                       MoreArgs = list(sl_folds = sl_training_folds,
                       all_fits = all_fits, 
                       learners = learners),
                       SIMPLIFY = FALSE)
    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_out, Y_out, pred_out)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred"); o})

    return(out)
}

#' Create input list for \code{get_risk_sl}
#' 
#' @param split_Y The outcome matrix split by relevant validation folds.
#' @param Yname The names of the outcomes. Used to search \code{all_fits}
#' and \code{all_sl}.
#' @param all_fits List of all learner fits.
#' @param all_fit_tasks List of all learner fit tasks (faster to search over than
#' search over all_fits). 
#' @param all_sl List of all super learner weight fits. 
#' @param all_weight List of all outcome weight fits. 
#' @param V Number of folds.
#' @param learners Vector of super learner vectors.
#' @param sl_control List of super learner controls.
#' @param folds Cross-validation folds.
#' 
#' @return List with each entry a list with entries: valid_fold (the number of the
#' corresponding fold), Y (univariate outcome in this validation fold), pred (matrix
#' of outcome predictions for this validation fold from learners fit in 
#' training folds). \code{get_risk_sl} is only used to compute the cross-validated 
#' risk of individual outcomes and so is only used in the outer most cross-validation
#' layer. 

get_risk_sl_input <- function(split_Y, Yname, all_fits, learners, V,
                              folds, all_sl, all_fit_tasks, all_weight, sl_control){
    # we need to get sl fit for this guy
    train_matrix <- utils::combn(V, V-1)

    pred_list <- lapply(split(train_matrix, col(train_matrix)), get_sl_pred_out, 
                    Ynames = Yname, outer_valid_folds = NULL, V = V,
                    all_fit_tasks = all_fit_tasks, all_fits = all_fits, 
                    all_sl = all_sl, sl_control = sl_control,
                    learners = learners)

    Y_list <- lapply(split(train_matrix, col(train_matrix)), get_Y_out, V = V, 
                split_Y = split_Y, training_folds = NULL)

    fold_list <- lapply(split(train_matrix, col(train_matrix)), 
                   get_fold_out, split_Y = split_Y, training_folds = NULL,
                   valid_folds = NULL, V = V)

    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_list, Y_list, pred_list)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred"); o})

    return(out)
}
#' Create input list for \code{get_risk_learner}
#' 
#' @param split_Y The outcome matrix split by relevant validation folds.
#' @param Yname The names of the outcomes. Used to search \code{all_fits}
#' and \code{all_sl}.
#' @param all_fits List of all learner fits.
#' @param all_fit_tasks List of all learner fit tasks (faster to search over than
#' search over all_fits). 
#' @param all_sl List of all super learner weight fits. 
#' @param all_weight List of all outcome weight fits. 
#' @param V Number of folds.
#' @param learner Name of super learner wrapper.
#' @param sl_control List of super learner controls.
#' @param folds Cross-validation folds.
#' 
#' @return List with each entry a list with entries: valid_fold (the number of the
#' corresponding fold), Y (univariate outcome in this validation fold), pred (matrix
#' of outcome predictions for this validation fold from learners fit in 
#' training folds). \code{get_risk_learner} is only used to compute the cross-validated 
#' risk of individual outcomes and so is only used in the outer most cross-validation
#' layer. 

get_risk_learner_input <- function(split_Y, Yname, learner, all_fits, V,
                              folds, all_sl, all_fit_tasks, all_weight, sl_control){
    # we need to get sl fit for this guy
    train_matrix <- utils::combn(V, V-1)

    pred_list <- lapply(split(train_matrix, col(train_matrix)), get_learner_pred_out, 
                    Ynames = Yname, outer_valid_folds = NULL, V = V,
                    all_fit_tasks = all_fit_tasks, all_fits = all_fits, 
                    sl_control = sl_control,
                    learner = learner)

    Y_list <- lapply(split(train_matrix, col(train_matrix)), get_Y_out, V = V, 
                split_Y = split_Y, training_folds = NULL)

    fold_list <- lapply(split(train_matrix, col(train_matrix)), 
                   get_fold_out, split_Y = split_Y, training_folds = NULL,
                   valid_folds = NULL, V = V)

    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_list, Y_list, pred_list)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred"); o})

    return(out)
}

#' Create input list for \code{get_risk}
#' 
#' \code{get_risk} computes the cross-validated risk of the entire 
#' procedure by calling \code{y_weight_control$cv_risk_fn} with this
#' input list. 
#' 
#' @param split_Y The outcome matrix split by relevant validation folds.
#' @param Ynames The names of the outcomes. Used to search \code{all_fits}
#' and \code{all_sl}.
#' @param all_fits List of all learner fits.
#' @param all_fit_tasks List of all learner fit tasks (faster to search over than
#' search over all_fits). 
#' @param all_sl List of all super learner weight fits. 
#' @param all_weight List of all outcome weight fits. 
#' @param V Number of folds.
#' @param learners Vector of super learner vectors.
#' @param sl_control List of super learner controls.
#' 
#' @return List with each entry a list with entries: valid_fold (the number of the
#' corresponding fold), Y (outcome matrix in this outer-most validation fold), 
#' pred (matrix of outcome predictions for this validation fold from super learner 
#' fit in V-1 training folds), and y_weight (vector of outcome weights computed by 
#' minimizing V-2 cross-validated risk of composite super learner). 
#' \code{get_risk_input} is only used to compute the cross-validated risk of 
#' composite super learner on the composite outcome 
#' and so is only used in the outer most cross-validation layer. 
#' 
# get V - 1 SL fits, get V - 1 SL predictions on validation, 
get_risk_input <- function(split_Y, Ynames, all_fits, 
                           V, learners, 
                           all_sl, all_fit_tasks, all_weight, sl_control){
    # we need to get sl fit for this guy
    train_matrix <- utils::combn(V, V-1)

    # need to search_fits_for_training_folds, but now for all Y
    # also need to search_sl_for_training_folds
    pred_list <- lapply(split(train_matrix, col(train_matrix)), get_sl_pred_out, 
                        Ynames = Ynames, outer_valid_folds = NULL, 
                        all_fit_tasks = all_fit_tasks, all_fits = all_fits, 
                        all_sl = all_sl, sl_control = sl_control,
                        V = V, learners = learners)

    # Y output
    Y_list <- apply(train_matrix, 2, get_Y_out, V = V, 
                    split_Y = split_Y, training_folds = NULL)
    
    # fold output 
    fold_list <- lapply(split(train_matrix, col(train_matrix)), 
                   get_fold_out, split_Y = split_Y, training_folds = NULL,
                   valid_folds = NULL, V = V)

    # outcome weights for each fold
    y_weight_list <- lapply(all_weight, "[", -1)

    # reorganize into proper format
    out <- do.call(Map, c(list, list(fold_list, Y_list, pred_list, y_weight_list)))
    out <- lapply(out, function(o){ names(o) <- c("valid_folds", "Y", "pred", "y_weight"); o})

    return(out)
}