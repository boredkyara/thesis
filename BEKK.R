# -----------------------------------------
# Author:
# Kyara Lucas 530510kl
# Optimized for compute engine
# The population is to the sample as the sample is to the bootstrap samples.
# -----------------------------------------

library(mvtnorm)
library(mgarchBEKK)
library(expm)
source("GMV.R")
set.seed(1234)

# Generate Multivariate BEKK(1,1,1) Garch Series
#
# Input:
# K               dimension of series
# T               length of series
# A, B, C         garch(1,1) parameters for volatilities
#
# Output:
# A list with the following components:
# epsi        (T x K) data matrix
# H           H_1, ... , H_T conditional covariance (K x K x T) matrix
multiBEKKGarch <- function(K, Tt, A, B, C){
  Tt = Tt + 50 # We discard the first 50 values
  
  ## Initialize empty list of matrices
  H = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
  for (t in 1:Tt){
    H[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize 3D matrices for Ht
  }
  
  epsi = matrix(data=0, nrow=Tt, ncol=K)
  temp =  solve(diag(K*K) - t(kronecker(A,A)) - t(kronecker(B,B))) %*% as.vector(crossprod(C))
  H[1,,] = matrix(temp,K,K)
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


# Calculate bootstrap CI's for BEKK(1,1,1)
# 
# Input:
# data          A 'data' list as outputted by multiBEKKGarch() function
# B             Number of bootstrap replicates 
#
# Output:
# A list with the following components:
# intervals     Confidence intervals for each parameter
# par_fit       Estimated parameters
pascualCI <- function(data, B){
  y = data$e
  H = data$H
  Tt = dim(data$e)[1]
  K = dim(data$e)[2]
  
  ### 1. Estimate (A,B,C) using QML
  sink("/dev/null")
  bekk.fit = suppressWarnings(BEKK(eps=y, order=c(1,1),
                                   params = c(0.8, 0.3, 0.8, 0.1, 0.1, 0.1, 0.1, 1, 0, 0, 1), # mgarchBEKK::BEKK
                                   method="L-BFGS-B", verbose=FALSE))
  sink()
  
  ## 1. Estimate weights using GMV
  par_fit = GMV(H[Tt,,])
  
  ## 2. Estimate the conditional var/covariances
  H_hat = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
  for (t in 1:Tt){
    H_hat[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize empty 3D matrices for Ht
  }
  
  C_hat = bekk.fit$est.params[[1]]
  A_hat = bekk.fit$est.params[[2]]
  B_hat = bekk.fit$est.params[[3]]
  
  temp = solve(diag(K*K) - t(kronecker(A_hat,A_hat)) - t(kronecker(B_hat,B_hat))) %*% as.vector(crossprod(C_hat))
  H_hat[1,,]  = matrix(data=temp, nrow=K, ncol=K)
  
  for(t in 2:Tt){
    H_hat[t,,] = crossprod(C_hat) + t(A_hat) %*% y[t-1,] %*% t(y[t-1,]) %*% A_hat + t(B_hat) %*% H_hat[t-1,,] %*% B_hat
  }
  
  ## 3. Compute residuals by the standardized observations
  z_hat = matrix(data=0, nrow=Tt, ncol=K)
  for(t in 1:Tt){
    z_hat[t,]= y[t,] %*% solve(chol(H_hat[t,,])) 
  }
  
  ### Replicate B bootstrappies (in parallel)
  replicates <- matrix(0, nrow=B, ncol=K)
  results <- list()
  n_cores = (detectCores() - 1)
  
  ## START MCLAPPLY ##
  results <- mclapply(X=rep(1,B), mc.cores=n_cores, FUN=function(b){
    
    r_boot = matrix(data=0, nrow=Tt, ncol=K)
    for(t in 1:Tt){
      r_boot[t,] = chol(H_hat[t,,]) %*% z_hat[sample(nrow(z_hat),1),]
    }
    
    ## 4. Estimate coefficients on bootstrap replicates
    fit.boot = suppressWarnings(BEKK(eps=r_boot, order=c(1,1),
                                     params = c(0.8, 0.3, 0.8, 0.1, 0.1, 0.1, 0.1, 1, 0, 0, 1), # mgarchBEKK::BEKK
                                     method="L-BFGS-B", verbose=FALSE))
    
    C_boot = fit.boot$est.params[[1]]
    A_boot = fit.boot$est.params[[2]]
    B_boot = fit.boot$est.params[[3]]
    
    # 5. Obtain boostrap replicates
    H_boot = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
    for (t in 1:Tt){
      H_boot[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize empty 3D matrices for H_boot
    }
    
    for(t in 1:Tt){
      if(t==1){
        H_boot[t,,] = H_hat[1,,]
      }else{
        H_boot[t,,] = crossprod(C_boot) + t(A_boot) %*% r_boot[t-1,] %*% t(r_boot[t-1,]) %*% A_boot + t(B_boot) %*% H_boot[t-1,,] %*% B_boot
      }
    }
    
    fit.boot.gmv = GMV(H=H_boot[Tt,,])
    fit.boot.gmv
  })
  ## END MCLAPPLY ##
  
  replicates = matrix(unlist(results), ncol=K, byrow=TRUE)
  
  # Confidence intervals
  intervals = matrix(data=0, nrow=K, ncol=2)
  for(i in 1:dim(replicates)[2]){
    intervals[i,] = c(unname(quantile(replicates[,i], c(0.05, 0.95))))
  }
  
  return(list("intervals"=intervals,
              "par_fit"=par_fit))
}

# Run simulations for coverage probability
# 
# Input:
# S               Number of simulations         
# Boot            Number of bootstrap replicates
# K, Tt, A, B, C  BEKK(1,1,1) DGP parameters
# trace           Boolean for print statements
#
# Output:
# A vector containing the coverage probabilites for each parameter
runSimulations <- function(S=10, Boot=10, K=2, Tt=100, A, B, C, trace=TRUE){
  nr_params = K
  
  start.time <- Sys.time()
  cat(sprintf("%s: starting simulation procedure ...\n", start.time))
  logFile = paste0("simulation ", start.time, ".txt")
  cat(paste0("Simulation started at ", start.time), file=logFile, append=FALSE, sep = "\n")
  cat("Coverage probability after every itteration: ", file=logFile, append=TRUE, sep = "\n")
  cat("---------------------------------------------", file=logFile, append=TRUE, sep = "\n")
  
  ## Initiate coverage count
  coverages = rep(0, nr_params)
  
  ## Simulation loop
  for(s in 1:S){
    ## Simulate BEKK(1,1,1) data
    data = multiBEKKGarch(K=K, Tt=Tt, A=A, B=B, C=C)
    
    ## Estimate CI's
    pascual_results = pascualCI(data=data, B=Boot)
    
    ## Update coverage count
    for(i in 1:nr_params){
      if(pascual_results$intervals[i,1] <= pascual_results$par_fit[i] && 
         pascual_results$intervals[i,2] >= pascual_results$par_fit[i]){
        coverages[i] = coverages[i] + 1
      }
    }
    
    cat(paste0(coverages/s), file=logFile, append=TRUE, sep = "\n")
    cat("---------------------------------------------", file=logFile, append=TRUE, sep = "\n")
    if(trace==TRUE){
      cat(sprintf("%s: %fth simulation completed.\n", Sys.time(), s))
    }
    
  }
  return(coverages/S)
}

## Set simulation DGP parameters
Boot = 999
K = 2
A = c(0.3, 0, -0.1, 0.3)
A = matrix(A, K, K)
B = c(0.9, -0.2, 0.2, 0.9)
B = matrix(B, K, K)
C = c(1, 0, 0.5, 1)
C = matrix(C, K, K)

test = runSimulations(S=20, Boot=999, K=2, Tt=800, A=A, B=B, C=C)


