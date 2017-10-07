context("Test Nonparametric R^2")

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

X_2 <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
Y_2 <- data.frame(Y1 = Y1, Y2 = Y2)

test_that("Nonparametric R^2 with default settings works", {
  
  fit <- cvma(Y = Y, X = X, V = 10, learners = c("SL.glm","SL.mean"))
  expect_equal(fit$cv_assoc$cv_measure, 0.768355, tolerance = 0.01)
})

test_that("Nonparametric R^2 with AUC SL risk works", {
  
  fit <- cvma(Y = Y_2, X = X_2, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              sl_control = list(ensemble_fn = "ensemble_linear",
                                optim_risk_fn = "optim_risk_sl_auc",
                                weight_fn = "weight_sl_01",
                                cv_risk_fn = "cv_risk_sl_r2",
                                family = gaussian(),
                                alpha = 0.05))
  
  expect_equal(fit$cv_assoc$cv_measure, -0.05150894, tolerance = 0.01)
})




