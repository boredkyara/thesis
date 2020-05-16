# -----------------------------------------
# Author:
# Kyara Lucas 530510kl
# The population is to the sample as the sample is to the bootstrap samples.
# -----------------------------------------

library(quadprog)
set.seed(1234)


# Calculate GMV (Global Minimum Variance) Portfolio Weights
#
# Input:
# H               Covariance matrix
#
# Output:
# Portfolio weights as vector
GMV <- function(H){
  nA = dim(H)[1]
  aMat  = array(1, dim = c(1,nA))
  bVec  = 1
  zeros = array(0, dim = c(nA,1))
  solQP <- solve.QP(H, zeros, t(aMat), bVec, meq = 1)
  return(solQP$solution)
}
