#' Cross-validated mean squared-error for computing super learner
#' 
#' 
#' In general, the function passed to \code{sl_control$optim_risk} should expect a named list
#' of outcomes (Y) and predictions (pred) in validation folds and should return a criteria by
#' which super learner weights may be optimized. The weights are input to the function via
#' \code{sl_weight} and are optimized in the \code{sl_control$weight_fn}. See
#' Examples section below for an example of the format of the \code{input} list used
#' for \code{sl_control$optim_risk} functions. 
#' 
#' In this case, the function computes cross-validated mean squared-error. 
#' 
#' @param sl_weight A numeric vector of super learner weights corresponding to each
#' \code{learner}. Typically, this is what is maximized over in \code{sl_control$weight_fn}.
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control Super learner control options. 
#' @export
#' 
#' @return Numeric value of cross-validated mean squared-error
#' @examples
#' 
#' # simulate data with proper format
#' # Y is one component of the multivariate outcome
#' # pred is the predictions made by learners 
#' input <- list(Y = rnorm(50), pred = cbind(rnorm(50), rnorm(50)))
#' 
#' # made up weights
#' sl_weight <- c(0.5, 0.5)
#' 
#' # linear ensemble
#' sl_control <- list(ensemble_fn = "ensemble_linear")
#' 
#' # get risk 
#' risk <- optim_risk_sl_se(sl_weight, input, sl_control)

optim_risk_sl_se <- function(sl_weight, input, sl_control){
    ens_pred <- do.call(sl_control$ensemble_fn, args = list(weight = sl_weight, pred = input$pred))
    risk <- mean((input$Y - ens_pred)^2)
    return(risk)
}

#' Cross-validated negative log-likelihood for computing super learner
#' 
#' In general, the function passed to \code{sl_control$optim_risk} should expect a named list
#' of outcomes (Y) and predictions (pred) in validation folds and should return a criteria by
#' which super learner weights may be optimized. The weights are input to the function via
#' \code{sl_weight} and are optimized in the \code{sl_control$weight_fn}. See
#' Examples section below for an example of the format of the \code{input} list used
#' for \code{sl_control$optim_risk} functions. 
#' 
#' In this case, the function computes cross-validated (quasi-) log-likelihood for the outcomes
#' scaled to be between 0 and 1. The option \code{trim} must be a value greater than
#' zero in order that the loss is bounded. The bounds on the outcome are set via
#' \code{l} and \code{u}. 
#' 
#' @param sl_weight A numeric vector of super learner weights corresponding to each
#' \code{learner}. Typically, this is what is maximized over in \code{sl_control$weight_fn}.
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control Super learner control options. 
#' @param l Lower bound on outcomes 
#' @param u Upper bound on outcomes
#' @param trim Where to trim \code{qlogis} to avoid \code{NA}s. 
#' 
#' @export
#' 
#' @examples
#' 
#' @return Numeric value of cross-validated negative log-likelihood
#' 
#' # simulate data with proper format
#' # Y is one component of the multivariate outcome
#' # pred is the predictions made by learners 
#' input <- list(Y = rnorm(50), pred = cbind(rnorm(50), rnorm(50)))
#' 
#' # made up weights
#' sl_weight <- c(0.5, 0.5)
#' # linear ensemble
#' sl_control <- list(ensemble_fn = "ensemble_linear")
#' # get risk setting l and u to min and max
#' # observed Y 
#' all_y <- unlist(lapply(input,"[",1))
#' risk <- optim_risk_sl_nloglik(sl_weight, input, l = min(all_y), u = max(all_y))
#' 

optim_risk_sl_nloglik <- function(sl_weight, input, sl_control, 
                             l = 0 , u = 1, trim = 0.001){
    ens_pred <- do.call(sl_control$ensemble_fn, 
                        args = list(weight = sl_weight, pred = input$pred))
    scale_ens_pred <- trim_p((ens_pred - l)/(u - l))
    scale_y <- trim_p((input$Y - l)/(u - l))
    risk <- - mean(scale_y * log(scale_ens_pred) + (1 - scale_y) * log(1 - scale_ens_pred))
    return(risk)
}



#' Cross-validated area under receiver operating characteristic curve (AUC)
#' 
#' 
#' In general, the function passed to \code{sl_control$optim_risk} should expect a named list
#' of outcomes (Y) and predictions (pred) in validation folds and should return a criteria by
#' which super learner weights may be optimized. The weights are input to the function via
#' \code{sl_weight} and are optimized in the \code{sl_control$weight_fn}. See
#' Examples section below for an example of the format of the \code{input} list used
#' for \code{sl_control$optim_risk} functions. 
#' 
#' In this case, the function computes cross-validated area under the receiver operating
#' characteristics curve (AUC) using \code{\link[cvAUC]{cvAUC}}. 
#' scaled to be between 0 and 1. The option \code{trim} must be a value greater than
#' zero in order that the loss is bounded. The bounds on the outcome are set via
#' \code{l} and \code{u}. 
#' 
#' @param sl_weight A numeric vector of super learner weights corresponding to each
#' \code{learner}. Typically, this is what is maximized over in \code{sl_control$weight_fn}.
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control Super learner control options. 
#' @export
#' @importFrom cvAUC cvAUC
#' @return Numeric value of cross-validated AUC. 
#' 
#' @examples
#' 
#' # simulate data with proper format
#' # Y is one component of the multivariate outcome
#' # pred is the predictions made by learners
#' input <- list(Y = rbinom(50, 1, 0.5), 
#'               pred = cbind(runif(50,0,1), runif(50,0,1)))
#' 
#' # made up weights
#' sl_weight <- c(0.5, 0.5)
#' # linear ensemble
#' sl_control <- list(ensemble_fn = "ensemble_linear")

#' # get risk 
#' all_y <- unlist(lapply(input,"[",1))
#' risk <- optim_risk_sl_auc(sl_weight, input, sl_control)
#' 
# TO DO: Test that this works
optim_risk_sl_auc <- function(sl_weight, input, sl_control){
    ens_pred <- lapply(input, function(i){
        do.call(sl_control$ensemble_fn, args = list(weight = sl_weight, pred = i$pred))
    })
    risk <- 1 - cvAUC::cvAUC(ens_pred, input$Y)
    return(risk)
}

#' Cross-validated non-parametric R-squared of the super learner
#' 
#' In general, the function passed to \code{sl_control$cv_risk} should expect a list
#' of outcomes and predictions in validation folds. The function should return a 
#' list with names cv_measure, ci_low, ci_high, and p_value. The output of this function 
#' is returned irrespective of the names of the list; however, the names are 
#' necessary for \code{print} methods to work properly.
#' 
#' In this case, the confidence intervals are computed on the scale of log(MSE/Var) and 
#' back-transformed to the R-squared scale. Here, MSE is the cross-validated mean 
#' squared-error of the super learner for predicting a univariate outcome (one of the
#' components of \code{Y}) and Var is the cross-validated marginal mean of this outcome. The p-value
#' is for the one-sided hypothesis test that cross-validated R-squared equals zero against
#' the alternative that it is greater than zero.
#' 
#' @param input List where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control List of super learner control options. 
#' 
#' @export
#' @importFrom stats qnorm
#' @return A list with named entries cv_measure, ci_low, ci_high, and p_value. The list will
#' be returned by \code{max_assoc} irrespective of the named entries; however, the \code{print}
#' methods will only work if the function returns the above names. 
#' @examples
#' 
#' # simulate data and properly format
#' set.seed(1234)
#' input <- list(list(Y = rnorm(50), pred = rnorm(50)),
#'               list(Y = rnorm(50), pred = rnorm(50)))
#' 
#' # get risk 
#' risk <- cv_risk_sl_r2(input, alpha = 0.05)

cv_risk_sl_r2 <- function(input, sl_control){
    # mse
    mse_list <- unlist(lapply(input, FUN = function(i){
        mean((i$Y-i$pred)^2)
    }), use.names = FALSE)
    cv_mse <- mean(mse_list)

    # var
    all_y <- unlist(lapply(input, "[", "Y"), use.names = FALSE)
    ybar <- mean(all_y)
    y_var <- mean((all_y - ybar)^2)

    # cross-validated r2
    cv_r2 <- 1 - cv_mse/y_var

    # influence functions
    all_pred <- unlist(lapply(input, "[", "pred"), use.names = FALSE)
    ic_mse <- (all_y - all_pred)^2 - cv_mse
    ic_var <- (all_y - ybar)^2 - y_var
        
    grad <- matrix(c(1/cv_mse, -1/y_var), ncol = 1)
    ic <- rbind(ic_mse, ic_var)

    se_1mlog_cv_r2 <- as.numeric(sqrt(tcrossprod(crossprod(grad, ic)))/length(ic_mse))
        
    ci_low <- 1 - exp(
        log(cv_mse/y_var) + stats::qnorm(1-(sl_control$alpha/2)) * se_1mlog_cv_r2
    )
    
    ci_high <- 1 - exp(
        log(cv_mse/y_var) - stats::qnorm(1-(sl_control$alpha/2)) * se_1mlog_cv_r2
    )
    p_value <- pnorm(log(cv_mse/y_var)/se_1mlog_cv_r2)
    
    return(list(cv_measure = cv_r2, ci_low = ci_low, ci_high = ci_high, p_value = p_value))
}

#' Cross-validated area under the receiver operating characteristic curve
#' for predictions based on super learner
#' 
#' In general, the function passed to \code{sl_control$cv_risk} should expect a list
#' of outcomes and predictions in validation folds. The function should return a 
#' list with names cv_measure, ci_low, ci_high, and p_value. The output of this function 
#' is returned irrespective of the names of the list; however, the names are 
#' necessary for \code{print} methods to work properly.
#' 
#' In this case, the confidence intervals are computed using the \code{\link[cvAUC]{cvAUC.ci}}
#' function from the \code{\link[cvAUC]{cvAUC}} package. The p-value
#' is for the one-sided hypothesis test that cross-validated AUC equals 0.5 against
#' the alternative that it is greater than 0.5.  
#' 
#' @param input List where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control List of super learner control options. 
#' 
cv_risk_sl_auc <- function(input, sl_control){
 
    # TO DO : Write cv_risk_sl_auc

}

#' Convex ensemble weights for super learner
#' 
#' In general, the function passed to \code{sl_control$weight_fn} should expect a list of 
#' named lists of outcomes (Y), predictions (pred) in validation folds. Typically,
#' this function is used to maximize \code{sl_control$optim_risk_fn} over 
#' weights. The function should return a named list. One of the names in the list should
#' be \code{weight}, which is the optimized weights. Other entries in the return list
#' are passed on to \code{sl_control$cv_risk_fn} (e.g., things needed to compute
#' cross-validated measure of association, though none are present for this particular 
#' function).
#'
#' In this case, the function uses \code{\link[Rsolnp]{solnp}} to minimize 
#' \code{sl_control$optim_risk_fn} over convex weights.
#'  
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control Super learner control options. 
#' 
#' @examples
#' 
#' # simulate data with proper format
#' input <- list(Y = rnorm(50), pred = cbind(rnorm(50), rnorm(50)))
#' # linear ensemble minimizing mean squared-error
#' sl_control <- list(ensemble_fn = "ensemble_linear",
#'                    optim_risk_fn = "optim_risk_sl_mse")
#' # get weights to minimize optim_risk
#' sl_weight <- weight_sl_convex(input, sl_control)
#' 
#' @export
#' 
#' @return Numeric vector giving convex weights for super learner. 

weight_sl_convex <- function(input, sl_control){
    # number of learners
    M <- dim(input[[1]]$pred)[2]

    # constrain weights to sum to 1
    constraint <- function(weight, input, sl_control){
        1 - sum(weight)
    }

    # reformat input into a single data.frame
    all_pred <- Reduce(rbind, lapply(input, '[[', "pred"))
    all_y <- Reduce(c, lapply(input, '[[', "Y"))
    solnp_input <- list(Y = all_y, pred = all_pred)

    # find weights
    fit <- Rsolnp::solnp(pars=rep(1/M, M), fun = eval(parse(text = sl_control$optim_risk_fn)), 
                         LB=rep(0, M), UB=rep(1, M), eqfun = constraint, 
                         control = list(trace=0), eqB=0, input = solnp_input, 
                         sl_control = sl_control)
    return(list(weight = fit$pars))
}

#' 0/1 weights for super learner (i.e., discrete super learner)
#' 
#' In general, the function passed to \code{sl_control$weight_fn} should expect a list of 
#' named lists of outcomes (Y), predictions (pred) in validation folds. Typically,
#' this function is used to maximize \code{sl_control$optim_risk_fn} over 
#' weights. The function should return a named list. One of the names in the list should
#' be \code{weight}, which is the optimized weights. Other entries in the return list
#' are passed on to \code{sl_control$cv_risk_fn} (e.g., things needed to compute
#' cross-validated measure of association, though none are present for this particular 
#' function).
#'
#' In this case, the function selects the single outcome with the lowest value 
#' returned by \code{sl_control$optim_risk_fn}
#' 
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (univariate outcome for this validation fold), pred (matrix of predictions
#' from \code{learner} and columns correspond to different \code{leaner}s). 
#' @param sl_control Super learner control options. 
#' @return Numeric vector giving 0/1 weights for super learner. 
#' @export
#' @examples
#' 
#' # simulate data and properly format
#' input <- list(Y = rnorm(50), pred = cbind(rnorm(50), rnorm(50)))
#' # linear ensemble
#' sl_control <- list(ensemble_fn = "ensemble_linear")
#' # get weights to minimize optim_risk
#' sl_weight <- weight_sl_01(input, sl_control)
#'  

weight_sl_01 <- function(input, sl_control){
    # number of learners
    M <- dim(input[[1]]$pred)[2]

    # reformat input into a single data.frame
    all_pred <- Reduce(rbind, lapply(input, '[[', "pred"))
    all_y <- Reduce(c, lapply(input, '[[', "Y"))
    risk_input <- list(Y = all_y, pred = all_pred)

    # get risk for each outcome
    risks <- rep(Inf, M)
    for(m in 1:M){
        y_weight <- rep(0, M)
        y_weight[m] <- 1
        risks[m] <- do.call(sl_control$optim_risk_fn, 
                            args = list(weight_y = weight_y, 
                                        input = risk_input, 
                                        sl_control = sl_control))    
    }

    # find lowest risk
    final_weight <- rep(0, M)
    final_weight[which.min(risks)] <- 1
    return(list(weight = final_weight))
}
