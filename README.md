<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`cvma`

[![Travis-CI Build Status](https://travis-ci.org/benkeser/cvma.svg?branch=master)](https://travis-ci.org/benkeser/cvma)
[![AppVeyor Build  Status](https://ci.appveyor.com/api/projects/status/github/benkeser/cvma?branch=master&svg=true)](https://ci.appveyor.com/project/benkeser/cvma)
[![Coverage Status](https://img.shields.io/codecov/c/github/benkeser/cvma/master.svg)](https://codecov.io/github/benkeser/cvma?branch=master)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Machine learning-based summary of association with multivariate outcomes

__Authors:__ [David Benkeser](https://www.benkeserstatistics.com/) and [Ivana Malenica](https://github.com/podTockom)

## Introduction

This package provides a method for summarizing the strength of association between a set of variables and a multivariate outcome. In particular, cross-validation is combined with stacked regression (aka super learning) to estimate the convex combination of a multivariate outcome that maximizes cross-validated R-squared of a super learner-based prediction. The method is particularly well suited for situations with high-dimensional covariates and/or complex relationships between covariates and outcomes. 

## Installation

You can install a stable release of `cvma` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:


```r
devtools::install_github("benkeser/cvma")
```

In the future, the package will be available from [CRAN](https://cran.r-project.org/) via


```r
install.packages("cvma")
```

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/benkeser/cvma/issues).

## Example

This minimal example shows how to use `cvma` with very simple, simulated data set. For more examples and detailed explanations, we refer the user to the vignette. To start with, we use the nonparametric R^2 to evaluate the strength of association between a set of variables and a multivariate outcome:




```r
suppressMessages(library(cvma))
set.seed(1234)

#Simulate data:
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rnorm(100, X$x1 + X$x2, 1)
Y2 <- rnorm(100, X$x1 + X$x2, 3)
Y <- data.frame(Y1 = Y1, Y2 = Y2)

#cvma with nonparametric R^2:
fit <- cvma(Y = Y, X = X, V = 10, 
                learners = c("SL.glm","SL.mean"))
fit
#>   cv_measure    ci_low   ci_high      p_value
#> 1  0.7648338 0.6616586 0.8365464 3.120828e-15
```

The following example evaluates the strength of association using AUC:


```r

#Simulate data:
X <- data.frame(x1=runif(n=100,0,5), x2=runif(n=100,0,5))
Y1 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1 + 0.2*X$x2))
Y2 <- rbinom(100, 1, plogis(-2 + 0.1*X$x1))
Y <- data.frame(Y1 = Y1, Y2 = Y2)

#cvma with AUC:
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
fit
#>   cv_measure    ci_low   ci_high   p_value
#> 1  0.3379673 0.1823169 0.4936178 0.9793412
```

## Variable importance

The cross-validated performance of two fits can be compared using the `compare_cvma` function. This can be used to define a variable importance measure for a set of variables. 


```r

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
# difference in cross-validated R^2 for the two fits
compare_cvma(fit1, fit2)
#>    contrast    ci_low   ci_high      p_value
#> 1 0.4587962 0.3195792 0.5980131 1.052884e-10
```

## License
&copy; 2017 [David C. Benkeser](http://www.benkeserstatistics.com)

The contents of this repository are distributed under the MIT license. See
below for details:
```
The MIT License (MIT)

Copyright (c) 2016-2017 David C. Benkeser

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
