context("Test negative log likelihood functions")

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

#Simulate some data:
set.seed(1234)

X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
Y <- data.frame(Y1 = Y1, Y2 = Y2)

test_that("Negative log-likelihood with default settings works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"),
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_nloglik",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_nloglik",
                                family = binomial(),
                                alpha = 0.05),
              y_weight_control = list(ensemble_fn = "ensemble_linear",
                                      weight_fn = "weight_y_01",
                                      optim_risk_fn = "optim_risk_y_nloglik",
                                      cv_risk_fn = "cv_risk_y_nloglik",
                                      alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.6731192, tolerance = 0.01)
})

test_that("Negative log-likelihood with mean squared-error SL risk works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"),
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_se",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_nloglik",
                                family = binomial(),
                                alpha = 0.05),
              y_weight_control = list(ensemble_fn = "ensemble_linear",
                                      weight_fn = "weight_y_01",
                                      optim_risk_fn = "optim_risk_y_nloglik",
                                      cv_risk_fn = "cv_risk_y_nloglik",
                                      alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.582109, tolerance = 0.01)
})

test_that("Negative log-likelihood with AUC SL risk and 01 weight works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"),
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_auc",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_nloglik",
                                family = binomial(),
                                alpha = 0.05),
              y_weight_control = list(ensemble_fn = "ensemble_linear",
                                      weight_fn = "weight_y_01",
                                      optim_risk_fn = "optim_risk_y_nloglik",
                                      cv_risk_fn = "cv_risk_y_nloglik",
                                      alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.5507803, tolerance = 0.01)
})

test_that("Negative log-likelihood with AUC y risk risk works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"),
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_se",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_nloglik",
                                family = binomial(),
                                alpha = 0.05),
              y_weight_control = list(ensemble_fn = "ensemble_linear",
                                      weight_fn = "weight_y_01",
                                      optim_risk_fn = "optim_risk_y_auc",
                                      cv_risk_fn = "cv_risk_y_nloglik",
                                      alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.4877479, tolerance = 0.01)
})


