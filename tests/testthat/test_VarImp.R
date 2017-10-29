context("Test Variable Importance")

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

#Simulate data:
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)

# fit data with full X
fit1 <- cvma(Y = Y, X = X, V = 10, 
             learners = c("SL.glm","SL.mean"))
# fit data with only x1
fit2 <- cvma(Y = Y, X = X[, -2, drop = FALSE], V = 10, 
             learners = c("SL.glm","SL.mean"))

test_that("Test Variable Importance with difference in cross-validated maximal association measures ", {

  fit<-compare_cvma(fit1, fit2, contrast = "diff")
  
  expect_equal(fit$contrast, 0.4874491, tolerance = 0.01)
})

test_that("Test Variable Importance with ratio of R^2", {

  fit<-compare_cvma(fit1, fit2, contrast = "ratio")
  
  expect_equal(fit$contrast, 2.757304, tolerance = 0.01)
  expect_equal(fit$ci_low, 1.31494, tolerance = 0.01)
  expect_equal(fit$ci_high, 4.199667, tolerance = 0.01)
})

test_that("Test Variable Importance with ratio of R^2, CI symmetric on log-scale", {

  fit<-compare_cvma(fit1, fit2, contrast = "logratio")
  
  expect_equal(fit$contrast, 2.757304, tolerance = 0.01)
  expect_equal(fit$ci_low, 1.510909, tolerance = 0.01)
  expect_equal(fit$ci_high, 5.031888, tolerance = 0.01)
})
