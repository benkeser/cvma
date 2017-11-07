context("Test saving fits option...")

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
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)

test_that("Summary refit function works", {
  
  fit <- cvma(Y = Y, X = X, V = 5, 
              learners = c("SL.glm","SL.mean"), 
              return_control = list(outer_weight = TRUE,
                                    outer_sl = TRUE,
                                    inner_sl = FALSE,
                                    all_y = TRUE,
                                    all_learner_assoc = TRUE,
                                    all_learner_fits = TRUE))
  
  expect_equal(fit$cv_assoc$cv_measure, 0.7623672, tolerance = 0.01)
  
  #now refit:
  re_fit <- reweight_cvma(fit, Y = Y, X = X, 
                          sl_control = list(ensemble_fn = "ensemble_linear",
                                            optim_risk_fn = "optim_risk_sl_se",
                                            weight_fn = "weight_sl_01",
                                            cv_risk_fn = "cv_risk_sl_r2",
                                            family = gaussian(),
                                            alpha = 0.05),
                          y_weight_control = list(ensemble_fn = "ensemble_linear",
                                                  weight_fn = "weight_y_01",
                                                  optim_risk_fn = "optim_risk_y_r2",
                                                  cv_risk_fn = "cv_risk_y_r2",
                                                  alpha = 0.05))
  
  expect_equal(re_fit$cv_assoc$cv_measure, 0.7541366, tolerance = 0.01)
  
})





