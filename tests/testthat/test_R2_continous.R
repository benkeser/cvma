context("Test Nonparametric R^2, continous outcome")

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
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)

test_that("Nonparametric R^2 with default settings works", {
  
  fit <- cvma(Y = Y, X = X, V = 10, learners = c("SL.glm","SL.mean"))
  expect_equal(fit$cv_assoc$cv_measure, 0.7648338, tolerance = 0.01)
})

test_that("Nonparametric R^2 with continous outcome, nloglik SL risk and convex sl weights works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_nloglik",
                                weight_fn = "weight_sl_convex",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = gaussian(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.5753772, tolerance = 0.01)
})

test_that("Nonparametric R^2 with continous outcome, mean squared-error SL risk and 01 sl weights works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_se",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = gaussian(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.7686046, tolerance = 0.01)
})

test_that("Nonparametric R^2 with more learners works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean", "SL.glmnet", "SL.nnet"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_se",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = gaussian(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.7600102, tolerance = 0.01)
})
