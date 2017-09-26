devtools::document("~/Dropbox/R/cvma")
install.packages("~/Dropbox/R/cvma/", repos = NULL, type = "source")
q("no")


# example with nonparametric R^2
set.seed(1234)
library(SuperLearner)
library(future)
library(cvma)
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)
fit <- cvma(Y = Y, X = X, V = 10, 
                learners = c("SL.glm","SL.mean"))

# example with AUC
set.seed(1234)
library(SuperLearner)
library(future)
library(cvma)

X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
Y <- data.frame(Y1 = Y1, Y2 = Y2)
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
