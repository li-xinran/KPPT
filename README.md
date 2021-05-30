## This is an R package for Kernel-based Partial Permutation Test


### Installation
library(devtools)

install_github("li-xinran/KPPT")

### library the package and other required packages
library(PPT)

library(mvtnorm); library(psych); library(quadprog); library(Matrix); library(lmtest)

### explanation for the main function
?partial.perm.test

### Example with one dimensional covariates
n = 200

X = matrix(runif(n, -1, 1), ncol = 1)

Y = as.vector( X + rnorm(n) * sqrt(0.1) )

Z = rbinom(n, 1, 0.5) + 1

ppt1 <- partial.perm.test(Y, X, Z, kernel.list = list(name = "Gaussian", alpha = 0.5, singular = 10^(-5)), perm.max = 500, sig.level = 0.05, alternative = "GP.pseudo", bandwidth.choice = "MLE.all")

ppt1$p.val

### Example with two dimensional covariates
n = 200

X = matrix(runif(2*n, -1, 1), ncol = 2)

Y = as.vector( 0.5 * X[,1] + 0.5 * X[,2] + rnorm(n) * sqrt(0.1) )

Z = rbinom(n, 1, 0.5) + 1

ppt2 <- partial.perm.test(Y, X, Z, kernel.list = list(name = "Gaussian", alpha = c(0.5, 0.5), singular = 10^(-5)), perm.max = 500, sig.level = 0.05, alternative = "GP.pseudo", bandwidth.choice = "MLE.all")

ppt2$p.val
