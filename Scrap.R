# simulating data from DCC-GARCH(1, 1) DGP
simulateDCC <- function(a, b, Q, a0, A, B, nobs, ncut=1000){
  nobs <- nobs + ncut             # ncut is the number of observations to be removed
  ndim <- nrow(Q)
  
  const <- (1-a-b)*Q
  
  z <- diag(0, nobs, ndim)
  DCC <- diag(0, nobs, ndim^2)
  vecQ2 <- diag(0, nobs, ndim^2)
  Q2 <- Q         # initial value of Q0
  zz <- Q
  h <- diag(0, nobs, ndim)
  eps <- diag(0, nobs, ndim)
  ht <- rep(0, ndim)
  et2 <- rep(0, ndim)
  for(i in 1:nobs){
    Q2 <- const + a*zz + b*Q2
    invdQ2 <- diag(1/sqrt(diag(Q2)))
    D2 <- invdQ2%*%Q2%*%invdQ2
    DCC[i, ] <- as.vector(D2)
    vecQ2[i, ] <- as.vector(Q2)
    
    # zt <- mvrnorm(1, rep(0, ndim), D2)     # standardized residual sdampling from N(0, R_t) whereR_t = D2 is DCC at t
    #############################################    
    cholD2 <- t(chol(D2))
    zt <- drop(cholD2%*%rnorm(ndim))
    #############################################    
    zz <- zt%o%zt
    z[i, ] <- zt
    
    ht <- a0 + A%*%et2 + B%*%ht           # conditional variance at t
    h[i, ] <- drop(ht)
    eps[i, ] <- zt*sqrt(ht)               # simulated observation at t
    et2 <- eps[i, ]^2
  }
  
  # A character vector/matrix for naming parameters
  name.id <- as.character(1:ndim)
  namev <- diag(0, ndim, ndim)
  for(i in 1:ndim){
    for(j in 1:ndim){
      namev[i, j] <- paste(name.id[i], name.id[j], sep="")
    }
  }
  # naming rows
  rownames(DCC) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(vecQ2) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(z) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(h) <- c(rep(0, ncut), 1:(nobs-ncut))
  rownames(eps) <- c(rep(0, ncut), 1:(nobs-ncut))
  # naming columns
  colnames(DCC) <- paste("R", namev, sep="")
  colnames(vecQ2) <- paste("R", namev, sep="")
  colnames(z) <- paste("Series", name.id, sep="")
  colnames(h) <- paste("Series", name.id, sep="")
  colnames(eps) <- paste("Series", name.id, sep="")
  
  list(DCC=zoo(DCC[-(1:ncut), ]), z=zoo(z[-(1:ncut), ]), Q=zoo(vecQ2[-(1:ncut), ]), h=zoo(h[-(1:ncut), ]), eps=zoo(eps[-(1:ncut), ]))
}




# function [ r, H, R ] = generateData( K, T, theta, para, burnIn )
# %% generating simulation data from DCC model
# %   Syntax:
#   %       [ r, H, R ] = generateData( K, T, theta, para, burnIn )
#   %   Input:
#     %       K = dimension of time series
#     %       T = length of time series
#     %       theta = garch(1,1) parameters for volatilities
#     %       para = [alpha beta] dynamic correlation parameters
#     %       burnIn = how many initial data to discard
#     %   Output:
#       %       r = T-by-K data matrix
#       %       H: H_1,...,H_T conditional covariance, K*K*T matrix.
#       %       R: R_1,...,R_T conditional correlation, K*K*T matrix.
#       
#       %%
#         alpha = para(1);
#         beta = para(2);
#         T0 = T+burnIn;
#         
#         % initial value 
#         Q = eye(K);
#         Rt = eye(K);
#         mu0 = zeros(1,K);
#         h = ones(1,K);
#         
#         %temp = rand(K,K);
#         %temp = (temp + temp')/2;
# %S = eye(K) + 0.1*temp;
# S = eye(K);
# 
# H = zeros(K,K,T0);
# R = ones(K,K,T0);
# 
# r = zeros(T0,K);
# for t = 1:T0
#     epsilon = mvnrnd(mu0,Rt);
#     rt = epsilon.*sqrt(h);
#     r(t,:) = rt;
#     Q = (1-alpha-beta)*S + alpha*(epsilon'*epsilon) + beta*Q; % Q_{t+1}
#         Rt = diag(1./sqrt(diag(Q)))*Q*diag(1./sqrt(diag(Q))); % R_{t+1}
#         h = theta(1) + theta(2)*rt.^2 + theta(3)*h;  % D_{t+1}
#         
#         R(:,:,t) = Rt;
#         H(:,:,t) = diag(sqrt(h))*Rt*diag(sqrt(h));
#         end
#         
#         r = r((burnIn+1):T0,:);
#         R = R(:,:,(burnIn+1):T0);
#         H = H(:,:,(burnIn+1):T0);
#         
#         end






# class BEKK_dgp:
#   def __init__(self,pars,N,T):
#   A = np.asmatrix(np.reshape(pars[0:N*N],(N,N)))
# B = np.asmatrix(np.reshape(pars[N*N:2*N*N],(N,N)))
# C = np.asmatrix(np.zeros((N,N)))
# C[np.triu_indices(N)] = pars[2*N*N:len(pars)]
# Ht = np.zeros((N,N,T))
# epsi = np.zeros((N,T))
# Ht[:,:,0] =C@C.T + np.identity(N)
# epsi[:,0] = np.random.normal(0,1,N)@np.asmatrix(np.linalg.inv(scipy.linalg.sqrtm(Ht[:,:,0])))
# for t in range(1,T):
#   eps = np.asmatrix(epsi[:,t-1]).T
# H = np.asmatrix(Ht[:,:,t-1])
# Ht[:,:,t] = C@C.T + A@eps@eps.T@A.T + B@H@B.T
# epsi[:,t] = np.random.normal(0,1,N)@np.asmatrix(np.linalg.inv(scipy.linalg.sqrtm(Ht[:,:,t])))
# self.epsi = epsi
# self.Ht = Ht


# sim = simulateBEKK(series.count=2, T=1000, 
#                    order=c(1,1), 
#                    params=c(1, 0.5, 1, 0.3, 0, -0.1, 0.3, 0.9, -0.2, 0.2, 0.9))
# series.1 = sim$eps[[1]]
# series.2 = sim$eps[[2]]
# bekk.sim.fit = BEKK(eps=cbind(series.1, series.2), order=c(1,1),
#                     method="L-BFGS-B", 
#                     params=c(1, 0.5, 1, 0.3, 0, -0.1, 0.3, 0.9, -0.2, 0.2, 0.9))
# bekk.sim.fit$est.params
# plot(sim$eps[[1]], type="l")
# plot(sim$eps[[2]], type="l")