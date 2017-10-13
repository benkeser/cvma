#' ci.cvAUC_withIC
#' 
#' This function is nearly verbatim \link[cvAUC]{ci.cvAUC} from the cvAUC package. 
#' The only difference is that it additionally returns estimated influence functions, 
#' which allows for variable importance measures to be computed later.
#' 
#' @param predictions A vector, matrix, list, or data frame containing the predictions.
#' @param labels A vector, matrix, list, or data frame containing the true class labels. Must have the 
#' same dimensions as \code{predictions}.
#' @param label.order The default ordering of the classes can be changed by supplying 
#' a vector containing the negative and the positive class label (negative label first, 
#' positive label second).
#' @param folds If specified, this must be a vector of fold ids equal in length to \code{predictions} 
#' and \code{labels}, or a list of length V (for V-fold cross-validation) of vectors of indexes for 
#' the observations contained in each fold. The \code{folds} argument must only be specified if 
#' the \code{predictions} and \code{labels} arguments are vectors.
#' @param confidence number between 0 and 1 that represents confidence level.
#' 
#' @importFrom cvAUC AUC cvAUC
#' @importFrom data.table data.table
#' 
#' @return A list containing the following named elements: 
#' \item{cvAUC}{Cross-validated area under the curve estimate.}
#' \item{se}{Standard error.}
#' \item{ci}{A vector of length two containing the upper and lower bounds for the confidence interval.}
#' \item{confidence}{A number between 0 and 1 representing the confidence.}
#' \item{ic}{A vector of the influence function evaluated at observations.}
#' @export

ci.cvAUC_withIC <- function(predictions, labels, label.ordering = NULL, 
                            folds = NULL, confidence = 0.95) {
  
  # Pre-process the input
  clean <- .process_input(predictions = predictions, labels = labels, 
                          label.ordering = label.ordering, folds = folds,
                          ids = NULL, confidence = confidence)
  
  predictions <- clean$predictions  # Length-V list of predicted values
  labels <- clean$labels  # Length-V list of true labels
  pos <- levels(labels[[1]])[[2]]  # Positive class label
  neg <- levels(labels[[1]])[[1]]  # Negative class label
  n_obs <- length(unlist(labels))  # Number of observations
  
  # Inverse probability weights across entire data set
  w1 <- 1/(sum(unlist(labels) == pos)/n_obs)  # Inverse weights for positive class
  w0 <- 1/(sum(unlist(labels) == neg)/n_obs)  # Inverse weights for negative class

  # This is required to cleanly get past R CMD CHECK
  # https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  pred <- label <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL  

  .IC <- function(fold_preds, fold_labels, pos, neg, w1, w0) {
      n_rows <- length(fold_labels)
      n_pos <- sum(fold_labels == pos)
      n_neg <- n_rows - n_pos
      auc <- cvAUC::AUC(fold_preds, fold_labels)
      DT <- data.table::data.table(pred = fold_preds, label = fold_labels)
      DT <- DT[order(pred, -xtfrm(label))]
      DT[, `:=`(fracNegLabelsWithSmallerPreds, cumsum(label == 
                                                        neg)/n_neg)]
      DT <- DT[order(-pred, label)]
      DT[, `:=`(fracPosLabelsWithLargerPreds, cumsum(label == 
                                                       pos)/n_pos)]
      DT[, `:=`(icVal, ifelse(label == pos, w1 * (fracNegLabelsWithSmallerPreds - 
                                                    auc), w0 * (fracPosLabelsWithLargerPreds - auc)))]
      return(DT$icVal)
  }

  icOut <- mapply(FUN = .IC, SIMPLIFY = FALSE, fold_preds = predictions, 
    fold_labels = labels, MoreArgs = list(pos = pos, neg = neg, w1 = w1, w0 = w0))
  # not back-sorted
  ic <- unlist(icOut)
  # Estimated variance
  sighat2 <- mean(unlist(lapply(icOut, var)))
  se <- sqrt(sighat2/n_obs)  
  cvauc <- cvAUC::cvAUC(predictions, labels)$cvAUC
  z <- qnorm(confidence + (1 - confidence)/2)
  ci_cvauc <- c(cvauc - (z * se), cvauc + (z * se))
  ci_cvauc[1] <- ifelse(ci_cvauc[1] < 0, 0, ci_cvauc[1])  #Truncate CI at [0,1]
  ci_cvauc[2] <- ifelse(ci_cvauc[2] > 1, 1, ci_cvauc[2]) 
  
  return(list(cvAUC = cvauc, se = se, ci = ci_cvauc, confidence = confidence, ic = ic))
}

#' Unexported function from cvAUC package
#' @param predictions A vector, matrix, list, or data frame containing the predictions.
#' @param labels A vector, matrix, list, or data frame containing the true class labels. Must have the 
#' same dimensions as \code{predictions}.
#' @param label.order The default ordering of the classes can be changed by supplying 
#' a vector containing the negative and the positive class label (negative label first, 
#' positive label second).
#' @param folds If specified, this must be a vector of fold ids equal in length to \code{predictions} 
#' and \code{labels}, or a list of length V (for V-fold cross-validation) of vectors of indexes for 
#' the observations contained in each fold. The \code{folds} argument must only be specified if 
#' the \code{predictions} and \code{labels} arguments are vectors.
#' @param ids Vector of ids
#' @param confidence confidence interval level
#' @importFrom ROCR prediction
.process_input <- function (predictions, labels, label.ordering = NULL, folds = NULL, 
    ids = NULL, confidence = NULL) 
{
    .vec_to_list <- function(idxs, vec) {
        return(vec[idxs])
    }
    if (!is.null(folds)) {
        if (class(predictions) == "list" | class(labels) == "list") {
            stop("If folds is specified, then predictions and labels must both be vectors.")
        }
        if (length(predictions) != length(labels)) {
            stop("predictions and labels must be equal length")
        }
        if (is.vector(folds) && !is.list(folds)) {
            if (length(folds) != length(labels)) {
                stop("folds vector must be the same length as the predictions/labels vectors.")
            }
            else {
                fids <- as.list(unique(folds))
                folds <- lapply(fids, function(fid, folds) {
                  which(folds == fid)
                }, folds)
            }
        }
        else if (!is.list(folds)) {
            stop("If specifying the folds argument, folds must be a list\n of vectors of indices that correspond to each CV fold or a vector of fold numbers\n the same size as the predictions/labels vectors.")
        }
        else if (length(unlist(folds)) != length(labels)) {
            stop("Number of observations in the folds argument does not equal number of predictions/labels.")
        }
        predictions <- sapply(folds, .vec_to_list, vec = predictions)
        labels <- sapply(folds, .vec_to_list, vec = labels)
        if (length(labels) > length(unlist(labels))) {
            stop("Number of folds cannot exceed the number of observations.")
        }
    }
    pred <- ROCR::prediction(predictions = predictions, labels = labels, 
        label.ordering = label.ordering)
    predictions <- pred@predictions
    labels <- pred@labels
    if (!is.null(ids)) {
        if (is.list(ids)) {
            if (length(unlist(ids)) != length(unlist(labels))) {
                stop("ids must contain same number of observations as predictions/labels.")
            }
        }
        else if (is.vector(ids)) {
            if (is.null(folds)) {
                ids <- list(ids)
            }
            else {
                ids <- sapply(folds, .vec_to_list, vec = ids)
            }
        }
        else if (is.matrix(ids) | is.data.frame(ids)) {
            ids <- as.list(data.frame(ids))
        }
        else {
            stop("Format of ids is invalid.")
        }
        if (length(ids) > 1) {
            n_ids <- sum(sapply(ids, function(i) {
                length(unique(i))
            }))
            if (length(unique(unlist(ids))) != n_ids) {
                warning("Observations with the same id are currently spread across multiple folds.\nAll observations with the same id must be in the same fold to avoid bias.")
            }
        }
    }
    if (!is.null(confidence)) {
        if (is.numeric(confidence) && length(confidence) == 1) {
            if (confidence <= 0 | confidence >= 1) {
                stop("confidence value must fall within (0,1)")
            }
        }
    }
    return(list(predictions = predictions, labels = labels, folds = folds, 
        ids = ids))
}