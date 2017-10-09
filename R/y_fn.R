#' Cross-validated non-parametric R-squared for computing composite outcome weights
#' 
#' In general, the function passed to \code{y_weight_control$optim_risk} should expect a named list
#' of outcomes (Y) and predictions (pred) in validation folds and should return a criteria by
#' which outcome weights may be optimized. The weights are input to the function via
#' \code{y_weight} and are optimized in the \code{y_weight_control$weight_fn}. See
#' Examples section below for an example of the format of the \code{input} list used
#' for \code{y_weight_control$optim_risk} functions. 
#'  
#' In this case, the function computes cross-validated nonparametric R-squared. 
#' 
#' @param y_weight A numeric vector of weights corresponding to each
#' outcome. Typically, this is what is maximized over in \code{y_weight_control$weight_fn}.
#' @param input A list with named entries Y (matrix of outcomes for this validation fold) and 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' 
#' @return Numeric value of cross-validated R-squared
#' @examples
#' 

# # simulate data with proper format
# input <- list(list(Y = cbind(rnorm(50), rnorm(50)), 
#                    pred = cbind(rnorm(50), rnorm(50))))
# 
# # made up weights
# y_weight <- c(0.5, 0.5)
# 
# # linear combination of outcomes
# y_weight_control <- list(ensemble_fn = "ensemble_linear")
# 
# # get risk 
# risk <- optim_risk_y_r2(y_weight, input, y_weight_control)

#TO DO: Check this example

optim_risk_y_r2 <- function(y_weight, input, y_weight_control){
    # get ensemble of outcomes
    ens_y <- do.call(y_weight_control$ensemble_fn, args = list(weight = y_weight, pred = input$Y))
    # get ensemble of predictions
    ens_p <- do.call(y_weight_control$ensemble_fn, args = list(weight = y_weight, pred = input$pred))
    # get r2
    r2 <- 1 - mean((ens_y - ens_p)^2)/mean((ens_y - mean(ens_y))^2)
    risk <- 1 - r2
    return(risk)
}

#' Cross-validated nonparametric R-squared for evaluating maximum association. 
#' 
#' In general, the function passed to \code{y_weight_control$cv_risk} should expect a list
#' of outcomes and predictions in validation folds, in addition to a list called
#' \code{y_weight} that contains the outcome weights (computed in the training sample)
#' corresponding to this validation fold and any other information needed by 
#' \code{y_weight_control$cv_risk} (e.g., anything needed to compute confidence 
#' intervals -- in this case the marginal mean of the composite outcome in the 
#' training sample). The function should return a list with names cv_measure, ci_low,
#' ci_high, and p_value. The output of this function is returned irrespective of the
#' names of the list; however, the names are necessary for \code{print} methods to 
#' work properly.
#' 
#' In this case, the confidence intervals are computed on the scale of log(MSE/Var) and 
#' back-transformed to the R-squared scale. Here, MSE is the cross-validated mean 
#' squared-error of the composite super learner predicting the composite outcome 
#' and Var is the cross-validated marginal mean of the composite outcome. The p-value
#' is for the one-sided hypothesis test that cross-validated R-squared equals zero against
#' the alternative that it is greater than zero. 
#' 
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (matrix of outcomes for this validation fold), 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' @importFrom stats qnorm pnorm
#' @return List with named components cv_measure (cross-validated AUC), ci_low (lower
#' 100(1 - \code{y_weight_control$alpha})\% CI), ci_high (upper
#' 100(1 - \code{y_weight_control$alpha})\% CI), p_value (p-value of test of null
#' hypothesis that cvAUC = 0.5).
#' 
#' @examples
#' 

#TO DO: Check this example

# # simulate data with proper format
# input <- list(list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)), 
#                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#                    y_weight = list(weight = c(0.5, 0.5))),
#               list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)),
#                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#                    y_weight = list(weight = c(0.25, 0.75))))
# 
# # linear combination of outcomes
# y_weight_control <- list(ensemble_fn = "ensemble_linear")
# 
# # get risk 
# cv_risk <- cv_risk_y_r2(input, y_weight_control)

cv_risk_y_r2 <- function(input, y_weight_control){
    # get ensemble y
    ens_y <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = i$y_weight$weight, pred = i$Y))  
    })

    # get ensemble p
    ens_p <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = i$y_weight$weight, pred = i$pred))  
    })

    # get ensemble marginal mean 
    ens_ybar <- lapply(input, function(i){
        rep(i$y_weight$ybar, dim(i$pred)[1])
    })

    # mse
    mse_list <- unlist(mapply(y = ens_y, p = ens_p, FUN = function(y, p){
        mean((y-p)^2)
    }), use.names = FALSE)
    cv_mse <- mean(mse_list)

    # var
    var_list <- mapply(y = ens_y, ybar = ens_ybar, FUN = function(y, ybar){
        mean((y-ybar)^2)
    })
    cv_var <- mean(var_list)

    # cross-validated r2
    cv_r2 <- 1 - cv_mse/cv_var

    ic_mse <- (unlist(ens_y) - unlist(ens_p))^2 - cv_mse
    ic_var <- (unlist(ens_y) - unlist(ens_ybar))^2 - cv_var
        
    grad <- matrix(c(1/cv_mse, -1/cv_var), ncol = 1)
    ic <- rbind(ic_mse, ic_var)

    se_1mlog_cv_r2 <- as.numeric(sqrt(tcrossprod(crossprod(grad, ic)))/length(ic_mse))
        
    ci_low <- 1 - exp(
        log(cv_mse/cv_var) + stats::qnorm(1-(y_weight_control$alpha/2)) * se_1mlog_cv_r2
    )
    
    ci_high <- 1 - exp(
        log(cv_mse/cv_var) - stats::qnorm(1-(y_weight_control$alpha/2)) * se_1mlog_cv_r2
    )
    p_value <- stats::pnorm(log(cv_mse/cv_var)/se_1mlog_cv_r2)
    
    return(list(cv_measure = cv_r2, ci_low = ci_low, ci_high = ci_high, p_value = p_value))
}

#' Cross-validated area under the receiver operating characteristic curve 
#' for computing composite outcome weights
#' 
#' In general, the function passed to \code{y_weight_control$optim_risk} should expect a named list
#' of outcomes (Y) and predictions (pred) in validation folds and should return a criteria by
#' which outcome weights may be optimized. The weights are input to the function via
#' \code{y_weight} and are optimized in the \code{y_weight_control$weight_fn}. See
#' Examples section below for an example of the format of the \code{input} list used
#' for \code{y_weight_control$optim_risk} functions. 
#'  
#' @param y_weight A numeric vector of weights corresponding to each
#' outcome. Typically, this is what is maximized over in \code{y_weight_control$weight_fn}.
#' @param input A list with named entries Y (matrix of outcomes for this validation fold) and 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' @importFrom cvAUC cvAUC
#' 
#' @return Numeric value of cross-validated AUC.
#' 
#' @examples
#' 
#' #Simulate data with proper format:
#' input <- list(list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5), rbinom(50,1,0.5)), 
#' pred = cbind(runif(50,0,1), runif(50,0,1), runif(50,0,1))))
#' 
#' #Linear combination of outcomes:
#' y_weight_control <- list(ensemble_fn = "ensemble_linear")
#' 
#' #Example weights:
#' y_weight<-c(0,1,0)
#' 
#' #Get risk:
#' risk <- optim_risk_y_auc(y_weight, input, y_weight_control)

optim_risk_y_auc <- function(y_weight, input, y_weight_control){
    ens_y <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = y_weight, pred = i$Y))
    })
    ens_p <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = y_weight, pred = i$pred))
    })
    if(!all(unlist(ens_y) %in% c(0,1))){
        stop("risk_y_auc requires all composite outcome to be either 0 or 1")
    }
    auc <- cvAUC::cvAUC(ens_p, ens_y, input$folds)
    risk <- 1 - auc$cvAUC
    return(risk)
}


#' Cross-validated area under the receiver operating characteristic curve 
#' for evaluating maximum association. 
#' 
#' In general, the function passed to \code{y_weight_control$cv_risk} should expect a list
#' of outcomes and predictions in validation folds, in addition to a list called
#' \code{y_weight} that contains the outcome weights (computed in the training sample)
#' corresponding to this validation fold and any other information needed by 
#' \code{y_weight_control$cv_risk} (e.g., anything needed to compute confidence 
#' intervals). The function should return a list with names cv_measure, ci_low,
#' ci_high, and p_value. The output of this function is returned irrespective of the
#' names of the list; however, the names are necessary for \code{print} methods to 
#' work properly.
#' 
#' In this case, the confidence intervals are computed using the \code{cvAUC::cvAUC.ci}
#' function from the \code{cvAUC::cvAUC} package. The p-value
#' is for the one-sided hypothesis test that cross-validated AUC equals 0.5 against
#' the alternative that it is greater than 0.5.  
#' 
#' @param input A list where each entry corresponds to a validation fold. Each entry is a list
#' with entries: Y (matrix of outcomes for this validation fold), 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' @importFrom cvAUC ci.cvAUC
#' @importFrom stats pnorm
#' @return List with named components cv_measure (cross-validated AUC), ci_low (lower
#' 100(1 - \code{y_weight_control$alpha})\% CI), ci_high (upper
#' 100(1 - \code{y_weight_control$alpha})\% CI), p_value (p-value of test of null
#' hypothesis that cvAUC = 0.5).
#' 
#' @examples
#' 
 
# # simulate data with proper format
# input <- list(list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)), 
#                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#                    y_weight = list(weight = c(0.5, 0.5))),
#               list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)),
#                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#                    y_weight = list(weight = c(0.25, 0.75))))
# 
# # linear combination of outcomes
# y_weight_control <- list(ensemble_fn = "ensemble_linear")
# 
# # get risk 
# cv_risk <- cv_risk_y_auc(input, y_weight_control)

#TO DO: Check this example

cv_risk_y_auc <- function(input, y_weight_control){
    ens_y <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = i$y_weight$weight, pred = i$Y))
    })
    ens_p <- lapply(input, function(i){
        do.call(y_weight_control$ensemble_fn, args = list(weight = i$y_weight$weight, pred = i$pred))
    })
    if(!all(unlist(ens_y) %in% c(0,1))){
        stop("cv_risk_y_auc requires all composite outcome to be either 0 or 1")
    }
    cv_auc_fit <- cvAUC::ci.cvAUC(ens_p, ens_y, confidence = 1 - y_weight_control$alpha)
    # p-value of one-sided test that cvAUC = 0.5 vs. cvAUC > 0.5
    p_value <- stats::pnorm((cv_auc_fit$cvAUC - 0.5) / cv_auc_fit$se, lower.tail = FALSE)
    out <- list(cv_measure = cv_auc_fit$cvAUC, ci_low = cv_auc_fit$ci[1],
                ci_high = cv_auc_fit$ci[2], p_value = p_value)
    return(out)
}

#' Compute optimized convex weights for outcomes
#' 
#' In general, the function passed to \code{y_weight_control$weight_fn} should expect a list of 
#' named lists of outcomes (Y), predictions (pred) in validation folds. Typically,
#' this function is used to maximize \code{y_weight_control$optim_risk_fn} over 
#' weights. The function should return a named list. One of the names in the list should
#' be \code{weight}, which is the optimized weights. Other entries in the return list
#' are passed on to \code{y_weight_control$cv_risk_fn} (e.g., things needed to compute
#' cross-validated measure of association).
#' 
#' In this case, the function uses \code{Rsolnp::solnp} to minimize 
#' \code{y_weight_control$optim_risk_fn} 
#' over convex weights. In addition to the optimized weights, the function
#' returns the marginal mean of the composite outcome based on the optimized weights, 
#' which is used by \code{\link{cv_risk_y_r2}} to compute the cross-validated nonparametric
#' R-squared.  
#' 
#' @param input A list where each entry corresponds to a validation fold. Each entry of \code{input}
#' is a list with entries: Y (matrix of outcomes for this validation fold), 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' @importFrom Rsolnp solnp
#' @return List with named components weight (optimized weights) and ybar (marginal
#' mean outcome for optimized weights, used to compute the cross-validated nonparametric
#' R-squared). 
#' 
#' @examples
#' 
#' # simulate data with proper format
#' input <- list(list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)), 
#'                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#'                    y_weight = list(weight = c(0.5, 0.5))),
#'               list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)),
#'                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#'                    y_weight = list(weight = c(0.25, 0.75))))
#' 
#' # linear combination of outcomes
#' y_weight_control <- list(ensemble_fn = "ensemble_linear",
#'                          optim_risk_fn = "optim_risk_y_r2")
#' 
#' # get risk 
#' weight <- weight_y_convex(input, y_weight_control)

weight_y_convex <- function(input, y_weight_control){
    # number of outcomes
    J <- dim(input[[1]]$Y)[2]

    # constrain weights to sum to 1
    constraint <- function(y_weight, input, y_weight_control){
        1 - sum(y_weight)
    }

    # reformat input into a single data.frame
    all_pred <- Reduce(rbind, lapply(input, '[[', "pred"))
    all_y <- Reduce(rbind, lapply(input, '[[', "Y"))
    solnp_input <- list(Y = all_y, pred = all_pred)

    # find weights
    fit <- Rsolnp::solnp(pars=rep(1/J, J), fun = eval(parse(text=y_weight_control$optim_risk_fn)), 
                         LB=rep(0, J), UB=rep(1, J), eqfun = constraint, 
                         control = list(trace = 0), eqB = 0, input = solnp_input, 
                         y_weight_control = y_weight_control)

    # get ybar for the minimized weights
    final_ens_y <- do.call(y_weight_control$ensemble_fn, args = list(weight = fit$pars, pred = all_y))

    return(list(weight = fit$pars, ybar = mean(final_ens_y)))
}


#' Find single outcome with best risk
#' 
#' In general, the function passed to \code{y_weight_control$weight_fn} should expect a list of 
#' named lists of outcomes (Y), predictions (pred) in validation folds. Typically,
#' this function is used to maximize \code{y_weight_control$optim_risk_fn} over 
#' weights. The function should return a named list. One of the names in the list should
#' be \code{weight}, which is the optimized weights. Other entries in the return list
#' are passed on to \code{y_weight_control$cv_risk_fn} (e.g., things needed to compute
#' cross-validated measure of association).
#' 
#' In this case, the function searches over all outcomes for the single outcome that
#' minimizes \code{y_weight_control$optim_risk_fn}. 
#' 
#' @param input A list where each entry corresponds to a validation fold. Each entry of \code{input}
#' is a list with entries: Y (matrix of outcomes for this validation fold), 
#' pred (matrix of super learner predictions for each outcomes with columns corresponding to 
#' different outcomes). 
#' @param y_weight_control Composite outcome weight control options. 
#' @export
#' @importFrom Rsolnp solnp
#' @return List with named components weight (optimized weights) and ybar (marginal
#' mean outcome for optimized weights, used to compute the cross-validated nonparametric
#' R-squared). 
#' 
#' @examples
#' 
#' # simulate data with proper format
#' input <- list(list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)), 
#'                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#'                    y_weight = list(weight = c(0.5, 0.5))),
#'               list(Y = cbind(rbinom(50,1,0.5), rbinom(50,1,0.5)),
#'                    pred = cbind(runif(50,0,1), runif(50,0,1)),
#'                    y_weight = list(weight = c(0.25, 0.75))))
#' 
#' # linear combination of outcomes ok with 0/1 weights
#' y_weight_control <- list(ensemble_fn = "ensemble_linear",
#'                          optim_risk_fn = "optim_risk_y_auc")
#' 
#' # get risk 
#' weight <- weight_y_01(input, y_weight_control)

weight_y_01 <- function(input, y_weight_control){
    # number of outcomes
    J <- dim(input[[1]]$Y)[2]
    
    # reformat input into a single data.frame
    all_pred <- Reduce(rbind, lapply(input, '[[', "pred"))
    all_y <- Reduce(rbind, lapply(input, '[[', "Y"))
    all_folds <- unlist(lapply(input,function(l){ rep(l$valid_folds, length(l$Y[,1])) }))
    solnp_input <- list(Y = all_y, pred = all_pred)

    # get risk for each outcome
    risks <- rep(0, J)
    for(j in 1:J){
        y_weight <- rep(0, J)
        y_weight[j] <- 1
        risks[j] <- do.call(y_weight_control$optim_risk_fn, 
                            args = list(y_weight = y_weight, input = input,
                                        y_weight_control = y_weight_control))
    }

    # find lowest risk
    final_weight <- rep(0, J)
    final_weight[which.min(risks)] <- 1
    return(list(weight = final_weight))
}
