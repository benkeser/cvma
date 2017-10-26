#' Helper function to make a task list for computing the outer super learner
#' @keywords internal
make_outer_sl_task_list <- function(Ynames, V){
    all_training <- 1:V
    out <- sapply(Ynames, function(y){
        list(Yname = y, training_folds = 1:V)
    }, simplify = FALSE)
    names(out) <- NULL    
    return(out)
}

#' Helper function to make a task list for computing cv risks of learners
#' @keywords internal
make_outer_learner_task_list <- function(Ynames, V, learners){
    all_training <- 1:V
    out <- vector(mode = "list", length = length(Ynames)*length(learners))
    ct <- 0
    for(y in Ynames){
        for(l in learners){
            ct <- ct + 1
            out[[ct]] <- list(training_folds = all_training,
                              SL_wrap = l,
                              Yname = y)
        }
    }
    names(out) <- NULL    
    return(out)
}

#' Helper function to make a task list for computing outcome weights
#' @keywords internal 
make_y_weight_task_list <- function(V){
    out <- utils::combn(V, V-1, FUN = function(x){
        return(list(training_folds = x))
    }, simplify = FALSE)
    return(out)
}
#' Helper function to make a task list for computing super learners.
#' @keywords internal
make_sl_task_list <- function(Ynames, V, fold_fits = c(V-1, V-2)){
    # get super learner fit on all 1-way hold-outs and 2-way hold-outs
    if(length(Ynames) == 1){
        fold_fits <- V - 1
    }
    all_training <- NULL
    for(v in fold_fits){
        all_training <- c(all_training, utils::combn(V, v, simplify = FALSE))
    }
    full_param <- expand.grid(training_folds = all_training, 
                              Yname = Ynames)
    # reorganize into a list
    out <- apply(full_param, 1, function(x){
        list(training_folds = x$training_folds, Yname = x$Yname)
    })
    return(out)   
}
#' Helper function to make a task list for computing learner fits.
#' @keywords internal
make_fit_task_list <- function(Ynames, learners, V, return_outer_sl,
                               fold_fits = NULL){
    out <- list()
    if(is.null(fold_fits)){
        # if more than 1 learner, then assume superlearner in which case
        # we need V-3 fold learner fits
        if(length(learners) > 1 & length(Ynames) > 1){
            fold_fits <- c(V - 1, V - 2, V - 3)        
        }else if(length(learners) + length(Ynames) > 2){
            # there's either multiple outcomes (in which case CV is needed 
            # to determine outcome weights) or there's multiple learners (in which
            # case CV is needed to make super learners)
            fold_fits <- c(V - 1, V - 2)
        }else{
            # there's only one outcome and one learner
            fold_fits <- V - 1
        }
    }
    # if sl return is requested, we'll need learner fits on all data
    if(return_outer_sl){
        fold_fits <- c(V, fold_fits)
    }
    for(v in fold_fits){
        com <- utils::combn(V, v, simplify = FALSE)
        full_param <- expand.grid(SL_wrap = learners,
                                  com = com,
                                  Yname = Ynames)
        tmp <- apply(full_param, 1, function(x){
            a <- list(training_folds = x$com,
                      Yname = x$Yname,
                      SL_wrap = x$SL_wrap)
            return(a)
        })
        out <- c(out, tmp)
    }
    return(out)
}