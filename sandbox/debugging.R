library(cvma)
list_control_options()
library(SuperLearner)
library(future)
set.seed(1234)
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)
fit <- cvma(Y = Y, X = X, V = 5, 
                 learners = c("SL.glm","SL.mean","SL.gam"))
set.seed(1234)
fit3 <- cvma(Y = Y, X = X[,-2,drop=FALSE], V = 5, 
                learners = c("SL.glm"))

set.seed(1234)
fit2 <- cvma(Y = Y, X = X[,-2,drop=FALSE], V = 5, 
                learners = c("SL.glm","SL.mean"))
set.seed(1234)
fit4 <- cvma(Y = Y[,-2,drop=FALSE], X = X[,-2,drop=FALSE], V = 5, 
                learners = c("SL.glm","SL.mean"))
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# it seems that fit4 should be same as fit2$cv_assoc_all_y[[1]][1:4], 
# and the point estimates are the same, but for some reason the confidence 
# intervals are difference, which is strange... need to debug to understand
# why the discrepancy, as this should not be the case!!!!!!!!!!!!!!!!!!!!!!!!
# ACTUALLY -- seems that this is ok, because of how ybar is being computed
# in cv_r2_y vs. cv_r2_sl -- the latter is what should be used in univariate
# calls to cvma.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cvma:::compare_cvma(fit, fit2)
cvma:::compare_cvma(fit, fit2, contrast = "ratio")
cvma:::compare_cvma(fit, fit2, contrast = "logratio")

summary(fit, "weights")
summary(fit, "superlearner")
summary(fit, "outcomes")
summary(fit, "learners")

  library(cvma)
  set.seed(1234)

  #Simulate data:
  X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
  Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
  Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
  Y <- data.frame(Y1 = Y1, Y2 = Y2)

  #cvma with AUC:
  set.seed(1234)
  fit <- cvma(Y = Y, X = X, V = 5, 
                  learners = c("SL.glm","SL.mean"),
                  sl_control = list(ensemble_fn = "ensemble_linear",
                                     optim_risk_fn = "optim_risk_sl_nloglik",
                                     weight_fn = "weight_sl_convex",
                                     cv_risk_fn = "cv_risk_sl_auc",
                                     family = binomial(),
                                     alpha = 0.05),
                  y_weight_control = list(ensemble_fn = "ensemble_linear",
                                    weight_fn = "weight_y_01",
                                    optim_risk_fn = "optim_risk_y_auc",
                                    cv_risk_fn = "cv_risk_y_auc",
                                    alpha = 0.05))
#cvma with AUC:
set.seed(1234)
fit2 <- cvma(Y = Y[,-2,drop = FALSE], X = X, V = 5, 
                learners = c("SL.glm","SL.mean"),
                sl_control = list(ensemble_fn = "ensemble_linear",
                                   optim_risk_fn = "optim_risk_sl_nloglik",
                                   weight_fn = "weight_sl_convex",
                                   cv_risk_fn = "cv_risk_sl_auc",
                                   family = binomial(),
                                   alpha = 0.05),
                y_weight_control = list(ensemble_fn = "ensemble_linear",
                                  weight_fn = "weight_y_01",
                                  optim_risk_fn = "optim_risk_y_auc",
                                  cv_risk_fn = "cv_risk_y_auc",
                                  alpha = 0.05))
cvma:::compare_cvma(fit, fit2)
cvma:::compare_cvma(fit, fit2, contrast = "ratio")
cvma:::compare_cvma(fit, fit2, contrast = "logratio")


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# try out the saving fits option
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
set.seed(1234)
library(SuperLearner)
library(future)
library(cvma)
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)
# results for super learner and R^2 for convex 
# combination of outcomes
fit <- cvma(Y = Y, X = X, V = 5, 
            learners = c("SL.glm","SL.mean"), 
        return_control = list(outer_weight = TRUE,
                                  outer_sl = TRUE,
                                  all_y = TRUE,
                                  all_learner_assoc = TRUE,
                                  all_learner_fits = TRUE))

# re-weight with discrete super learner and R^2 
# for a single outcome
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
cvma:::compare_cvma(fit, re_fit)


#!!!!!!!!!!!!!!!!!!!!!!!
# checking that negative log likelihood functions work


  library(cvma)
  set.seed(1234)

  #Simulate data:
  X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
  Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
  Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
  Y <- data.frame(Y1 = Y1, Y2 = Y2)

#cvma with AUC:
set.seed(1234)
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
                              alpha = 0.05),
            return_control = list(outer_weight = TRUE,
                                  outer_sl = TRUE,
                                  inner_sl = TRUE, 
                                  all_y = TRUE,
                                  all_learner_assoc = TRUE,
                                  all_learner_fits = FALSE))
