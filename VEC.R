library(mvtnorm)
library(rmgarch)
library(MTS)
library(expm)

set.seed(1234)

## TODO:
# Generate VEC DGP
# Fit VEC model and check parameters

# Generate Multivariate VEC Garch Series
#
# Input:
# K               dimension of series
# T               length of series
# theta           garch(1,1) parameters for volatilities (omega, alpha, beta)
# alpha, beta     dynamic correlation parameters
#
# Output:
# r       (T x K) data matrix
# H       H_1, ... , H_T conditional covariance (K x K x T) matrix
# R       R_1, ... , R_T conditional correlation (K x K x T) matrix
multi_vec_garch <- function(K=2, Tt=100, a, b, c){
  A = a
  B = b
  C = c
  
  ## Initialize empty list of matrices
  H = array(rep(0, Tt*K*K), dim=c(Tt, K, K))
  for (t in 1:Tt){
    H[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize 3D matrices for Ht
  }
  
  Ht[1,,] = diag(K)
  vecht = array(rep(1, K*(K+1)/2*Tt), dim=c(K*(K+1)/2, Tt)) # vecht[,t] = c(1, 1, 1)
  epsi = matrix(data=0, nrow=K, ncol=Tt) # epsi[,t] = c(1, 1)
  
  ## Recursively generate VEC series
  for (t in 2:Tt){
    eta = vech(crossprod(epsi[,t-1]))
    ht = C + A%*%eta + B%*%H[t-1,,]
    
  }
  
  
}

# for:
#   d = C + A@vec(np.asmatrix(epsi[:,t-1])@np.asmatrix(epsi[:,t-1]).T)+ B@vec(np.asmatrix(Ht[:,:,t - 1]))
#   vecht[:,t] = d.A1
#   Ht[:,:,t] = restack(np.asmatrix(vecht[:,t]).T)
#   epsi[:,t-1] = scipy.linalg.sqrtm(Ht[:,:,t])@np.random.normal(0,1,(N,1))
# return np.asmatrix(epsi)