# -----------------------------------------
# Author:
# Kyara Lucas 530510kl
# -----------------------------------------

## TODO:
# Generate BEKK DGP x
# Fit BEKK model and check parameters x

## Import libraries
library(mvtnorm)
library(mgarchBEKK)
library(expm)
set.seed(1234)

# Generate Multivariate BEKK(1,1,K) Garch Series
#
# Input:
# K               dimension of series
# T               length of series
# A, B, C         garch(1,1) parameters for volatilities
#
# Output:
# A list with the following components:
# epsi       (T x K) data matrix
# H       H_1, ... , H_T conditional covariance (K x K x T) matrix
multiBEKKGarch <- function(K=2, Tt=100, A, B, C){
  Tt = Tt + 50 # We discard the first 50 values
  
  ## Initialize empty list of matrices
  H = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
  for (t in 1:Tt){
    H[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize 3D matrices for Ht
  }
  
  epsi = matrix(data=0, nrow=Tt, ncol=K)
  H[1,,] = (C %*% t(C)) + diag(K)
  epsi[1,] = chol(H[1,,]) %*% matrix(rmvnorm(n=1, mean=rep(0, K), sigma=diag(K)))
  
  ## Recursively generate series
  for(t in 2:Tt){
    H[t,,] = crossprod(C) + (t(A) %*% epsi[t-1,] %*% t(epsi[t-1,]) %*% A) + (t(B) %*% H[t-1,,] %*% B)
    epsi[t,] = chol(H[t,,]) %*% matrix(rmvnorm(n=1, mean=rep(0,K), sigma=diag(K)))
  }
  
  ## We discard the first 50 values
  return(list("e"=epsi[50:Tt,],
              "H"=H[50:Tt,,]))
}

## Simulate BEKK(1,1,K) series
# Set parameters
K = 2
Tt = 1000
A = c(0.3, 0, -0.1, 0.3)
A = matrix(A, K, K)
B = c(0.9, -0.2, 0.2, 0.9)
B = matrix(B, K, K)
C = c(1, 0, 0.5, 1)
C = matrix(C, K, K)

bekk.sim = multiBEKKGarch(K=K, Tt=Tt, A=A, B=B, C=C)
bekk.fit = suppressWarnings(BEKK(eps=bekk.sim$e, order=c(1,1),
                params = c(0.8, 0.3, 0.8, 0.1, 0.1, 0.1, 0.1, 1, 0, 0, 1), # mgarchBEKK::BEKK
                method="L-BFGS-B", verbose=FALSE))
bekk.fit$est.params

plot(bekk.sim$e[,1], type="l")
plot(bekk.sim$e[,2], type="l")

