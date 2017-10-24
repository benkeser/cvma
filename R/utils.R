#' Print the output of \code{cvma}. Only prints 
#' the cross-validated summary measures. To print more
#' see \code{?summary.cvma}.
#' 
#' @param x An object of class \code{cvma}
#' @param ... Other options (not currently used)
#' 
#' @export
print.cvma <- function(x, ...){
	print(data.frame(x$cv_assoc[c("cv_measure","ci_low","ci_high","p_value")]))
}


#' Summarize results of \code{cvma} fit. Can be used to summarize
#' the weights obtained for each outcome (in both outer and inner
#' layers of cross-validation), the super learner fit for each outcome
#' (with learners fit on all the data), or the cross-validated associations
#' with each outcome. 
#' 
#' @param object An object of class \code{cvma}
#' @param aspect Can be \code{"weights"}, \code{"superlearner"},
#' \code{"outcomes"}, or \code{"learners"}
#' @param ... Other options (not currently used)
#' 
#' @export
#' 
#' @return A list summarizing relevant aspects of the \code{cvma} object
#' 
summary.cvma <- function(object, aspect = "outcomes", ...){
	if(aspect == "outcomes"){
		if(is.null(object$cv_assoc_all_y)){
			stop("To summarize associations with each outcome, return_all_y must be TRUE in call to cvma.")
		}
		out <- lapply(object$cv_assoc_all_y, function(cva){
			tmp <- data.frame(cva[c("cv_measure","ci_low","ci_high","p_value")])
		})
		names(out) <- object$y_names
		return(out)
	}else if(aspect == "weights"){
		out <- list()
		if(!is.null(object$outer_weight)){
			out$outer_weights <- data.frame(matrix(NA, nrow = 1,
		                                       ncol = length(object$y_names) + 1))
			colnames(out$outer_weights) <- c("training_folds", object$y_names)
			out$outer_weights$training_folds <- list(object$outer_weight$training_folds)
			out$outer_weights[,2:ncol(out$outer_weights)] <- object$outer_weight$weight
		}
		out$inner_weights <- data.frame(matrix(NA, nrow = length(object$inner_weight),
		                                       ncol = length(object$y_names) + 1))
		colnames(out$inner_weights) <- c("training_folds", object$y_names)
		out$inner_weights$training_folds <- lapply(object$inner_weight, "[[", "training_folds")
		out$inner_weights[,2:ncol(out$inner_weights)] <- Reduce(rbind, lapply(object$inner_weight, "[[", "weight"))
		return(out)
	}else if(aspect == "superlearner"){
		if(is.null(object$cv_assoc_all_y)){
			stop("To summarize super learner, return_outer_sl must be TRUE in call to cvma.")
		}
		out <- lapply(object$sl_fits, function(fit){
			tmp <- data.frame(fit[c("learner_names","learner_risks","sl_weight")])
		})
		names(out) <- object$y_names
		return(out)
	}else if(aspect == "learners"){
		if(is.null(object$cv_assoc_all_learners)){
			stop("To summarize learners, return_all_learners must be TRUE in call to cvma.")
		}
		n_learners <- length(object$cv_assoc_all_learners) / length(object$y_names)
		tmp_out <- split(object$cv_assoc_all_learners, sort(rep(1:length(object$y_names), n_learners)))
		names(tmp_out) <- object$y_names
		out <- lapply(tmp_out, function(tmpo){
			tmp <- lapply(tmpo, function(cva){
				data.frame(cva[c("SL_wrap","cv_measure","ci_low","ci_high","p_value")])
			})
			tmp2 <- Reduce(rbind, tmp)
			tmp2 <- tmp2[order(-tmp2$cv_measure),]
			row.names(tmp2) <- NULL
			tmp2
		})
		return(out)
	}
}

