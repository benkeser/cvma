context("Test summary functons...")

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

test_that("Summary functions are as expected...", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean","SL.gam"))
  
  expect_equal(summary(fit, "weights")[[1]]$Y1, 0.8801622, tolerance = 0.01)
  expect_equal(summary(fit, "weights")[[1]]$Y2, 0.1198378, tolerance = 0.01)
  expect_equal(fit$outer_weight$weight, c(0.8801622,0.1198378), tolerance = 0.01)
  expect_equal(summary(fit, "outcomes")[[1]]$cv_measure, 0.7334427, tolerance = 0.01)
  expect_equal(summary(fit, "learners")[[1]]$cv_measure[1], 0.7366201, tolerance = 0.01)
  expect_equal(summary(fit, "superlearner")[[1]]$learner_risks[1], 1.13373, tolerance = 0.01)
  
})


