#' Compare measures of association between two fits
#' 
#' @param object1 An object of class \code{cvma}
#' @param object2 An object of class \code{cvma}
#' @param contrast A character indicating what contrast to take. Can be \code{"diff"},
#' \code{"ratio"}, or \code{"logratio"}.
#' @param compare_control Contains information regarding format of the contrast we are interested in, its inverse function, 
#' type of contrast, IC components of the corresponding covariance calculation, and null hypothesis. 
#' Defaults to functions necessary for \code{"diff"}.
#' @param alpha Confidence intervals of 1 - alpha. 
#' @export
#' @importFrom stats cov
#' @return Estimates, CI, and p-value.
#' @examples
#' set.seed(1234)
#' library(SuperLearner)
#' library(future)
#' X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
#' Y1 <- rnorm(100, X$x1 + X$x2, 1)
#' Y2 <- rnorm(100, X$x1 + X$x2, 3)
#' Y <- data.frame(Y1 = Y1, Y2 = Y2)
#' fit1 <- cvma(Y = Y, X = X, V = 5, 
#'                 learners = c("SL.glm","SL.mean"))
#' fit2 <- cvma(Y = Y, X = X[,-2,drop = FALSE], V = 5, 
#'                 learners = c("SL.glm","SL.mean"))
#' # difference in R^2
#' compare_cvma(fit1, fit2, "diff")
#' # ratio of R^2, symmetric CI
#' compare_cvma(fit1, fit2, "ratio")
#' # ratio of R^2, CI symmetric on log-scale
#' compare_cvma(fit1, fit2, "logratio")

compare_cvma <- function(object1, object2, 
                         contrast = "diff",
                         compare_control = list(f = function(h_cv_measures){ h_cv_measures },
                                                f_inv = function(h_cv_measures){ h_cv_measures },
                                                h = function(cv_measure1, cv_measure2){
                                                	cv_measure1 - cv_measure2
                                                },
                                                fh_grad = function(cv_measure1, cv_measure2){
                                                	c(1, -1)
                                                }, null = 0),
                         alpha = 0.05){

	if(!is.null(contrast)){
		compare_control <- do.call(paste0("compare_control_",contrast), args = list())
	}else{
		if(!all(c("f","f_inv","h","fh_grad","null") %in% names(compare_control))){
			stop("some function missing in contrast. see ?compare_cvma for help.")
		}
	}

	this_c <- do.call(compare_control$h, args=list(cv_measure1 = object1$cv_assoc$cv_measure,
	                                       cv_measure2 = object2$cv_assoc$cv_measure))
	f_this_c <- do.call(compare_control$f, args=list(h_cv_measures = this_c))

	grad <- matrix(do.call(compare_control$fh_grad,
	                       args=list(cv_measure1 = object1$cv_assoc$cv_measure,
	                                 cv_measure2 = object2$cv_assoc$cv_measure)),nrow = 2)
	
	v <- stats::cov(cbind(object1$cv_assoc$ic, object2$cv_assoc$ic))
	n <- length(object1$cv_assoc$ic)

	se <- sqrt(t(grad)%*%v%*%grad / n)
	transform_ci <- rep(f_this_c, 3) + stats::qnorm(c(0.5,alpha/2,1- alpha/2))*
		rep(se,3)
	ci <- do.call(compare_control$f_inv, args = list(h_cv_measures = transform_ci))
	p_value <- 2*pnorm(-abs((f_this_c - compare_control$null)/se))
	
	return(cbind.data.frame(contrast = ci[1], ci_low = ci[2], ci_high = ci[3], p_value = p_value))
	#return(list(contrast = ci[1], ci_low = ci[2], ci_high = ci[3], p_value = p_value))
}

#' compare_control for taking a difference in cross-validated maximal 
#' association measures.
compare_control_diff <- function(){
	list(f = function(h_cv_measures){ h_cv_measures },
	     f_inv = function(h_cv_measures){ h_cv_measures },
	     h = function(cv_measure1, cv_measure2){
	     	cv_measure1 - cv_measure2
	     },
	     fh_grad = function(cv_measure1, cv_measure2){
	     	c(1, -1)
	     }, 
	     null = 0)
}

#' compare_control for taking a ratio of cross-validated maximal 
#' association measures.
compare_control_ratio <- function(){
	list(f = function(h_cv_measures){ h_cv_measures },
	     f_inv = function(h_cv_measures){ h_cv_measures },
	     h = function(cv_measure1, cv_measure2){
	     	cv_measure1/cv_measure2
	     },
	     fh_grad = function(cv_measure1, cv_measure2){
	     	c(1/cv_measure2, -cv_measure1/cv_measure2^2)
	     },
	     null = 1)
}

#' compare_control for taking a ratio of cross-validated maximal 
#' association measures, but computing CI and hypothesis test on log scale. 
compare_control_logratio <- function(){
	list(f = function(h_cv_measures){ log(h_cv_measures) },
	     f_inv = function(h_cv_measures){ exp(h_cv_measures) },
	     h = function(cv_measure1, cv_measure2){
	     	cv_measure1/cv_measure2
	     },
	     fh_grad = function(cv_measure1, cv_measure2){
	     	c(1/cv_measure2, -1/cv_measure2)
	     },
	     null = 0)
}
