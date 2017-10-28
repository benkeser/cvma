context("Test Nonparametric R^2, binary outcome")

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

X_2 <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1_2 <- rbinom(100, 1, plogis(-2 + 0.1*X_2$x1 + 0.2*X_2$x2))
Y2_2 <- rbinom(100, 1, plogis(-2 + 0.1*X_2$x1))
Y_2 <- data.frame(Y1 = Y1_2, Y2 = Y2_2)

test_that("Nonparametric R^2 with AUC SL risk and 01 sl weight works", {
  
  fit <- cvma(Y = Y_2, X = X_2, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_auc",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = binomial(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, -0.003636513, tolerance = 0.01)
})

test_that("Nonparametric R^2 with nloglik SL risk and convex sl weight works", {
  
  fit <- cvma(Y = Y_2, X = X_2, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_nloglik",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = binomial(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.001180872, tolerance = 0.01)
})


test_that("Nonparametric R^2 with nloglik SL risk and 01 sl weights works", {
  
  fit <- cvma(Y = Y_2, X = X_2, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_nloglik",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = binomial(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, -0.01554361, tolerance = 0.01)
})


test_that("Nonparametric R^2 with 01 y weight and AUC risk works ", {
  
  fit <- cvma(Y = Y_2, X = X_2, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_se",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = gaussian(),
                                alpha = 0.05),
              y_weight_control = list(ensemble_fn = "ensemble_linear",
                                      weight_fn = "weight_y_01",
                                      optim_risk_fn = "optim_risk_y_auc",
                                      cv_risk_fn = "cv_risk_y_r2",
                                      alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.0248, tolerance = 0.01)
})
