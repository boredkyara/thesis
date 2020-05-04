library(mvtnorm)
library(rmgarch)
library(expm)

set.seed(1234)

## TODO:
# - Generate DCC-Garch(1,1) series x
# - Estimate parameters DCC-Garch(1,1) series (check if matches) x
# - Implement bootstrap estimation
# - Implement CI and coverage probabilities

# Generate Multivariate DCC Garch Series
#
# Input:
# K               dimension of series
# T               length of series
# theta           garch(1,1) parameters for volatilities (omega, alpha, beta)
# alpha, beta     dynamic correlation parameters

# Output:
# r       (T x K) data matrix
# H       H_1, ... , H_T conditional covariance (K x K x T) matrix
# R       R_1, ... , R_T conditional correlation (K x K x T) matrix
multi_dcc_garch <- function(K=2, Tt=100, theta=c(0.015, 0.94), alpha=0.1, beta=0.85, omega=0.1){
  burn_in = Tt/5
  Tt = (6/5) * Tt
  
  a = theta[1] # Dynamic correlation alpha
  b = theta[2] # Dynamic correlation beta
  
  ## Set initial values
  Q = diag(K)
  Rt = diag(K)
  mu0 = matrix(data=0, nrow=1, ncol=K)
  h = matrix(data=1, nrow=1, ncol=K)
  
  Q_bar = diag(K) # cov(e, e')
  
  ## Initialize empty list of matrices
  H = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
  for (t in 1:Tt){
    H[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize 3D matrices for Ht
  }
  
  R = array(rep(1, Tt*K*K), dim=c(Tt, K, K))
  for(t in 1:Tt){
    R[t,,] = matrix(data=1, nrow=K, ncol=K) # Initialize 3D matrices for Rt
  }
  
  r = matrix(data=1, nrow=Tt, ncol=K)
  e = matrix(data=1, nrow=Tt, ncol=K)
  
  ## Generate DCC-GARCH(1,1) process recursively
  for(t in 1:Tt){
    epsilon = rmvnorm(1, mean=mu0, sigma=Rt)
    rt = epsilon*sqrt(h)
    
    r[t,] = rt
    e[t,] = epsilon
    
    # Q = (1-alpha-beta)*S + alpha*(epsilon'*epsilon) + beta*Q; % Q_{t+1}
    Q = (1-a-b)*Q_bar + a*(crossprod(epsilon)) + b*Q # Q_t+1
    
    # Rt = diag(1./sqrt(diag(Q)))*Q*diag(1./sqrt(diag(Q))); % R_{t+1}
    Qinvd = diag(1/sqrt(diag(Q)))
    Rt = Qinvd %*% Q %*% Qinvd
    
    h = omega + alpha * rt^2 + beta * h # h = sigma squared 
    Dt = diag(c(sqrt(h))) # D_t+1
    
    R[t,,] = Rt
    H[t,,] = Dt %*% Rt %*% Dt
  }
  
  r = r[burn_in:Tt,]
  R = R[burn_in:Tt,,]
  H = H[burn_in:Tt,,]

  return(list("r"=r,
              "R"=R,
              "H"=H,
              "e"=e))
}

simulated_data = multi_dcc_garch()



specifications = multispec(replicate(2, ugarchspec(mean.model = list(armaOrder = c(0,0)))))
dcc = dccspec(uspec = specifications, dccOrder = c(1, 1), distribution = 'mvnorm')
dcc_fit = dccfit(dcc, data = simulated_data$r)
dcc_fit
