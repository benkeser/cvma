#' Compute a linear ensemble of predictions
#' 
#' In general, \code{sl_control$ensemble_fn} and \code{y_weight_control$ensemble_fn}
#' must accept named inputs, pred and weight, which are, respectively, a matrix and a vector. 
#' The matrix has columns corresponding to different learner fits (for \code{sl_control$ensemble_fn})
#' or to the different components of the multivariate outcome (for \code{y_weight_control$ensemble_fn}).
#' 
#' In this case, the function computes a linear combination. 
#' 
#' @param pred A matrix of predictions
#' @param weight A vector of ensemble weights (should have same length as number
#' of columns in \code{pred})
#' @return A vector of ensembled values with length equal to the number of 
#' rows of \code{pred}. 
#' @export
ensemble_linear <- function(pred, weight){
    nonzero_idx <- weight != 0 
    weight <- matrix(weight, nrow = 1)[, nonzero_idx, drop = FALSE]
    out <- tcrossprod(weight, as.matrix(pred[ , nonzero_idx, drop = FALSE]))
    return(as.numeric(out))
}

#' Helper function to compute a trimmed logit
#' @keywords internal
#' @importFrom stats qlogis

trim_qlogis <- function(p, trim = 1e-4){
    return(stats::qlogis(trim_p(p, trim = trim)))
}

#' Helper function to trim a prediction
#' @keywords internal
trim_p <- function(p, trim = 1e-4){
    p[p < trim] <- trim; p[p > 1 - trim] <- 1 - trim
    return(p)
}


#' Compute a logit-linear ensemble of predictions
#' 
#' In general, \code{sl_control$ensemble_fn} and \code{y_weight_control$ensemble_fn}
#' must accept named inputs, pred and weight, which are, respectively, a matrix and a vector. 
#' The matrix has columns corresponding to different learner fits (for \code{sl_control$ensemble_fn})
#' or to the different components of the multivariate outcome (for \code{y_weight_control$ensemble_fn}).
#' 
#' In this case, the function computes an ensemble of scaled (by \code{l} and \code{u})
#' \code{pred} on the logit scale and back-transforms to the original scale. The option
#' \code{trim} is used to avoid numeric instabilities. 
#' 
#' @param pred A matrix of predictions
#' @param weight A vector of ensemble weights (should have same length as number
#' of columns in \code{pred})
#' @param l The lower scaling factor.
#' @param u The upper scaling factor
#' @param trim A value to trim logit computations to avoid numerical instabilities. 
#' @return A vector of ensembled values with length equal to the number of 
#' rows of \code{pred}. 
#' @export
#' @importFrom stats plogis
ensemble_logit_linear <- function(pred, weight, l = 0, u = 1, 
                             trim = 0.001){
    nonzero_idx <- weight != 0 
    scale_pred <- (pred - l)/(u - l)
    logit_scale_pred <- trim_qlogis(scale_pred, trim = trim)
    lin_pred <- tcrossprod(weight[nonzero_idx, ], logit_scale_pred[ , nonzero_idx, drop = FALSE])
    out <- stats::plogis(lin_pred)*(u - l) + l
    return(as.numeric(out))
}

