#' Helper function to make a task list for computing the outer super learner
#' @keyword internal
make_outer_sl_task_list <- function(Ynames, V){
    all_training <- 1:V
    out <- sapply(Ynames, function(y){
        list(Yname = y, training_folds = 1:V)
    }, simplify = FALSE)
    names(out) <- NULL    
    return(out)
}

#' Helper function to make a task list for computing outcome weights
#' @keyword internal 
make_y_weight_task_list <- function(V){
    out <- combn(V, V-1, FUN = function(x){
        return(list(training_folds = x))
    }, simplify = FALSE)
    return(out)
}
#' Helper function to make a task list for computing super learners
#' @keyword internal
make_sl_task_list <- function(Ynames, V){
    # get super learner fit on all 1-way hold-outs and 2-way hold-outs
    all_training <- c(combn(V, V-1, simplify = FALSE), combn(V, V-2, simplify = FALSE))
    full_param <- expand.grid(training_folds = all_training, 
                              Yname = Ynames)
    # reorganize into a list
    out <- apply(full_param, 1, function(x){
        list(training_folds = x$training_folds, Yname = x$Yname)
    })
    return(out)   
}
#' Helper function to make a task list for computing learner fits
#' @keyword internal
make_fit_task_list <- function(Ynames, learners, V, return_outer_sl){
    out <- list()
    fold_fits <- c(V-1,V-2,V-3)
    # if sl return is requested, we'll need learner fits on all data
    if(return_outer_sl){
        fold_fits <- c(V, fold_fits)
    }
    for(v in fold_fits){
        com <- combn(V, v, simplify = FALSE)
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