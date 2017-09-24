#' Helper function to search fits for particular outcomes and training_folds.
#' @keywords internal
search_fits_one_y <- function(l, Yname, training_folds){
    l$Yname == Yname & identical(l$training_folds, training_folds)
}
#' Helper function to search fits for particular training_folds. 
#' @keywords internal
search_fits_all_y <- function(l, training_folds){
    identical(l$training_folds, training_folds)
}
#' Helper function to search fits for particular outcomes (or single outcomes)
#' and training_folds.
#' @keywords internal
search_fits_for_training_folds <- function(training_folds, fits, y = NULL){
    # get the index of corresponding fits
    if(!is.null(y)){
        fit_idx <- which(sapply(fits, FUN = search_fits_one_y, 
                                Yname = y, training_folds = training_folds))
    }else{
        fit_idx <- which(sapply(fits, FUN = search_fits_all_y, 
                                training_folds = training_folds))
    }
    return(fit_idx)
}