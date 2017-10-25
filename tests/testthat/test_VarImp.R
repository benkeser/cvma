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

test_that("Test Variable Importance function", {

  # fit data with full X
  fit1 <- cvma(Y = Y, X = X, V = 10, 
               learners = c("SL.glm","SL.mean"))
  # fit data with only x1
  fit2 <- cvma(Y = Y, X = X[, -2, drop = FALSE], V = 10, 
               learners = c("SL.glm","SL.mean"))
  
  # difference in cross-validated R^2 for the two fits
  fit<-compare_cvma(fit1, fit2)
  
  expect_equal(fit$contrast, 0.4874491, tolerance = 0.01)
})
