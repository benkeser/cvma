context("Test AUC, binary outcome")

if (FALSE) {
  setwd("..")
  setwd("..")
  getwd()
  library("devtools")
  document()
  load_all("./")  # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # devtools::check() # runs full check
  setwd("..")
  install("cvma", build_vignettes = FALSE, dependencies = FALSE)  # INSTALL W/ devtools:
}

library(testthat)
library(SuperLearner)
library(future)
library(cvma)

set.seed(1234)

X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
Y <- data.frame(Y1 = Y1, Y2 = Y2)

#test_that("AUC with weight 01 works", {
#  
#  fit <- cvma(Y = Y, X = X, V = 5, scale=FALSE,
#              learners = c("SL.glm","SL.mean"),
#              sl_control = list(ensemble_fn = "ensemble_linear",
#                                optim_risk_fn = "optim_risk_sl_se",
#                                weight_fn = "weight_sl_convex",
#                                cv_risk_fn = "cv_risk_sl_r2",
#                                family = gaussian(),
#                                alpha = 0.05),
#              y_weight_control = list(ensemble_fn = "ensemble_linear",
#                                      weight_fn = "weight_y_01",
#                                      optim_risk_fn = "optim_risk_y_auc",
#                                     cv_risk_fn = "cv_risk_y_auc",
#                                      alpha = 0.05))
#  
#  expect_equal(fit$cv_assoc$cv_measure, 0.7102452, tolerance = 0.01)
#})

test_that("AUC with weight convex throws an error", {
  
  sl_control = list(ensemble_fn = "ensemble_linear",optim_risk_fn = "optim_risk_sl_nloglik",weight_fn = "weight_sl_convex",cv_risk_fn = "cv_risk_sl_auc",family = binomial(),alpha = 0.05)
  y_weight_control = list(ensemble_fn = "ensemble_linear",weight_fn = "weight_y_convex",optim_risk_fn = "optim_risk_y_auc",cv_risk_fn = "cv_risk_y_auc",alpha = 0.05)
  learners = c("SL.glm","SL.mean")
  
  expect_error(broken_fit <- cvma(Y,X,V=5,learners=learners,sl_control=sl_control,y_weight_control=y_weight_control),
               "risk_y_auc requires all composite outcome to be either 0 or 1.")
})

