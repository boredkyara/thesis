# -----------------------------------------
# Author:
# Kyara Lucas 530510kl
# -----------------------------------------

# Import libraries
library(fGarch)
set.seed(1234)

# Generate univariate GARCH(1,1) data
# 
# Input:
# Tnr           number of observations to be generated
# mu            mean parameter
# w             omega parameter
# a             alpha parameter
# b             beta parameter
# 
# Output:
# A list with the following components:
# y             series of generated y values
# sigma         series of generated sigma values
# e             series of generated residual values
# par_true      list of true parameter values

uniGarch <- function(Tnr, mu, w, a, b){
  par_true <- list("mu"=mu,"omega"=w, "alpha"=a, "beta"=b) # List of our model paramters
  sigma_t <- sqrt(w/(1-a-b)) # Set initial volatility
  
  sigma_seq <- c() # Initiate empty array for sequence of volatility
  y_seq <- c() # Initiate empty array for sequence of returns
  e_seq <- c() # Initiate empty array for sequence of white noise
  
  ## Generate GARCH(1,1) process recursively
  for(i in 1:Tnr){
    e_seq[i] = rnorm(1,0,1) # Generate white noise process w. unit variance
    # y_t1 = sigma_t * e_seq[i]
    y_t1 = rnorm(1,0,sigma_t)
    y_seq[i] = y_t1
    sigma_t1 = w + a*(y_t1^2) + b*(sigma_t^2)
    sigma_seq[i] = sigma_t1
    sigma_t = sqrt(sigma_t1)
  }
  
  return(list("y"=y_seq,
              "sigma"=sigma_seq,
              "e"=e_seq,
              "par_true"=par_true))
}


# Calculate bootstrap CI's for GARCH(1,1) as proposed in Pascual et. al (2006)
# 
# Input:
# data          A 'data' list as outputted by uni_garch() function
# B             Number of bootstrap replicates 
#
# Output:
# A list with the following components:
# intervals     Confidence intervals for each GARCH(1,1) parameter
# par_fit       Initial estimated parameter (using QML)
pascualCI <- function(data, B){
  
  ## 1. Estimate (w,a,b) using QML
  g = garchFit(formula = ~ garch(1,1), data = data$y, cond.dist="QMLE", trace=FALSE)
  par_fit = g@fit$coef
  par_fit = list("mu"=unname(par_fit[1]), "omega"=unname(par_fit[2]), 
                 "alpha"=unname(par_fit[3]), "beta"=unname(par_fit[4]))
  
  ## 2. Estimate conditional variances
  sigma2_seq_hat = c()
  sigma2_seq_hat[1] = (par_fit$omega/(1 - par_fit$alpha - par_fit$beta))
  
  Tnr = length(data$y)
  for(t in 2:Tnr){
    sigma2_hat = par_fit$omega + par_fit$alpha*data$y[t-1]^2 + par_fit$beta*sigma2_seq_hat[t-1]
    sigma2_seq_hat[t] = sigma2_hat
  }
  
  names(sigma2_hat) <- NULL
  # plot(sqrt(sigma2_seq_hat), type = "l")
  
  ## 3. Compute residuals
  eps_hat = data$y/sqrt(sigma2_seq_hat)
  # eps_hat = eps_hat - mean(eps_hat)
  
  replicates <- matrix(nrow = B, ncol = length(par_fit))
  
  # Replicate B bootstrappies
  for(b in 1:B){
    
    ## 4. Obtain bootstrap replicates {y_1^*, ..., y_T^*}
    sigma2_seq_boot = c()
    y_seq_boot = c()
    
    for(t in 1:Tnr){
      if(t==1){
        sigma2_boot = sigma2_seq_hat[1]
      } else{
        sigma2_boot = unname(par_fit$omega + par_fit$alpha*y_seq_boot[t-1]^2 
                             + par_fit$beta*sigma2_seq_boot[t-1]) 
      }
      y_boot = sqrt(sigma2_boot)*sample(eps_hat,1)
      
      sigma2_seq_boot[t] = sigma2_boot
      y_seq_boot[t] = y_boot
    }
    
    ## 5. Estimate coefficients on bootstrap replicates
    fit_boot = garchFit(formula = ~ garch(1,1), data = y_seq_boot, cond.dist="QMLE", trace=FALSE)
    replicates[b,] <- unname(fit_boot@fit$coef)
    
  }
  
  # Confidence intervals
  intervals = list("mu"=c(unname(quantile(replicates[,1], c(0.05, 0.95)))),
                   "omega"=c(unname(quantile(replicates[,2], c(0.05, 0.95)))),
                   "alpha"=c(unname(quantile(replicates[,3], c(0.05, 0.95)))),
                   "beta"=c(unname(quantile(replicates[,4], c(0.05, 0.95)))))
  
  return(list("intervals"=intervals,
              "par_fit"=par_fit))
}



runSimulation <- function(S, B, trace=TRUE){
  
  start.time <- Sys.time()
  cat(sprintf("%s: starting simulation procedure ...\n", start.time))
  
  # Initiate coverage count
  mu_coverage = 0
  omega_coverage = 0
  alpha_coverage = 0
  beta_coverage = 0
  
  for(s in 1:S){
    # Simulate univariate GARCH(1,1) data
    data = uniGarch(Tnr=500, mu=0, w=0.05, a=0.1, b=0.85)
    
    # Estimate CI's according to Pascual et. al (2006)
    pascual_results = pascualCI(data=data, B=B)
    
    # Update coverage count
    ## TODO: Change to fitted values
    if(pascual_results$intervals$mu[1] < pascual_results$par_fit$mu && 
       pascual_results$intervals$mu[2] > pascual_results$par_fit$mu){
      mu_coverage = mu_coverage + 1
    }
    
    if(pascual_results$intervals$omega[1] < pascual_results$par_fit$omega && 
       pascual_results$intervals$omega[2] > pascual_results$par_fit$omega){
      omega_coverage = omega_coverage + 1
    }
    
    if(pascual_results$intervals$alpha[1] < pascual_results$par_fit$alpha && 
       pascual_results$intervals$alpha[2] > pascual_results$par_fit$alpha){
      alpha_coverage = alpha_coverage + 1
    }
    
    if(pascual_results$intervals$beta[1] < pascual_results$par_fit$beta && 
       pascual_results$intervals$beta[2] > pascual_results$par_fit$beta){
      beta_coverage = beta_coverage + 1
    }
    
    if((s %% 10 == 0) && (trace==TRUE)){
      cat(sprintf("%s: %fth simulation completed.\n", Sys.time(), s))
    }
  }
  
  end.time <- Sys.time()
  cat(sprintf("%s: finished.\n", end.time))
  time.taken <- (end.time - start.time)
  cat(sprintf("%s: total minutes taken.\n", time.taken))
  
  return(list("mu_coverage"=mu_coverage/S,
              "beta_coverage"=beta_coverage/S,
              "omega_coverage"=omega_coverage/S,
              "alpha_coverage"=alpha_coverage/S))
}

results = runSimulation(S=100, B=500)
results

