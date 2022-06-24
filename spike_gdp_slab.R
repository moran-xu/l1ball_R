require(statmod)
require(truncnorm)

spike_gdp_slab<- function(y,X, steps =1000,
                          alpha = 0.5, eta_tilde = 1,
                          a_sigma=2,b_sigma= 1, 
                          a_w = 1,b_w = p,
                          burn =500){
  
  
  
  softThresholding<- function(beta,mu){
    sign(beta)*(abs(beta)-mu)*(abs(beta)>mu)
  }
  
  
  invSoftThresholding<- function(beta,theta,mu){
    S<- abs(theta)>0
    beta[S]<- sign(theta[S])*(abs(theta[S])+mu)
    beta
  }
  
  updateThetaS<- function(theta, S, sigma2,lam){
    
    X_s= as.matrix(X[,S])
    if(sum(S)>1){
      V = solve(t(X_s)%*%X_s + diag(1/tau[S]))
    }else{
      V= solve(t(X_s)%*%X_s + 1/tau[S])
    }
    mean_theta_s = V%*% ( t(X_s)%*%y)
    theta[S] <- t(chol(V))%*% rnorm(sum(S)) * sqrt(sigma2) + mean_theta_s
    theta
  }
  
  updateSigma2<- function(beta, tau, theta){
    
    S = sum(theta!=0)
    
    alpha_sigma<- (n+ sum(S))/2+a_sigma
    beta_sigma<- (sum(theta*theta/tau)+ sum((y-X%*%theta)**2))/2+ b_sigma
    sigma2 = 1/rgamma(1,alpha_sigma,rate = beta_sigma)
    if(sigma2>10){
      sigma2=10 # truncate at 10 to prevent numeric issue
    }
    sigma2
  }
  
  
  
  
  updateLambda<-function(theta,sigma2){
    rgamma( p, alpha+1, abs(theta)/sqrt(sigma2) + eta_tilde)
  }
  
  updateTau<- function(lam,sigma2,theta){
    1/rinvgauss(p, mean= sqrt(lam**2* sigma2/(theta**2)), shape =lam**2)
  }
  
  updateBetaNotS<- function(beta,tau,sigma2, mu){
    if(sum(!S)>0){
      beta[!S]<- runif(sum(!S), -mu,mu)
    }
    beta
  }
  
  updateW<- function(S){
    n_S = sum(S)
    rbeta(1, a_w+n_S, b_w+p-n_S)
  }
  
  logPost<- function(beta,sigma2,w){
    theta<- softThresholding(beta,mu)
    part1<-  - sum((y-X%*%theta)**2)/2/sigma2
    part2<- - (alpha+1) * sum(log(1+ abs(theta)/eta_tilde/sqrt(sigma2)))
    
    S = abs(beta)>mu
    
    part3 <- (a_w-1+ sum(S))*log(w)+(b_w-1+ p-sum(S))*log(1-w) 
    
    if(sum(abs(theta))==0){return(-Inf)}
    else{
      return(part1+ part2+ part3)
    }
  }
  
  
  randomWalkBeta<- function(beta,sigma2,w, eps=0.1,subset_size= 100){
    
    accept =0
    mh_count =0
    
    cur_beta<- beta
    curLogPost<- logPost(beta,sigma2,w)
    
    full_idx<- sample(c(1:p))
    
    
    idx0<- 1
    idx1<- idx0 + subset_size-1
    
    while(idx0<p){
      if(idx1>p)idx1=p
      beta[full_idx[idx0:idx1] ]<- beta[full_idx[idx0:idx1]]+ runif(idx1-idx0+1,-1,1)*eps
      propLogPost<- logPost(beta,sigma2,w)
      
      mh_count = mh_count+1
      
      if (log(runif(1))< propLogPost -curLogPost){
        curLogPost = propLogPost
        cur_beta[full_idx[idx0:idx1]] = beta[full_idx[idx0:idx1]]
        accept = accept + 1
      }else{
        beta[full_idx[idx0:idx1]] = cur_beta[full_idx[idx0:idx1]]
      }
      idx0 = idx0+ subset_size
      idx1 = idx1+ subset_size
    }
    return(list(beta, accept/mh_count))
  }
  
  
  n<- nrow(X)
  p<- ncol(X)
  
  beta<- runif(p)
  mu<- 0.5
  
  sigma2<- 0.1
  
  
  # theta<-  theta0 
  theta<- solve(t(X)%*%X + diag(0.1,p), t(X)%*%y)
  S<- theta!=0
  beta<- invSoftThresholding(beta,theta,mu)
  beta[!S]= mu/2
  lam <- rep(1,p)
  
  tau<- rep(1,p)
  
  
  eps_beta= 0.1
  w = 1/p
  
  trace_theta<- matrix(0,0,p)
  
  
  
  trace_theta<- matrix(0,0,p)
  trace_sigma2<- numeric(0)
  
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = steps, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for(step in c(1:steps)){
    lam<- updateLambda(theta,sigma2)
    tau<- updateTau(lam,sigma2,theta)
    theta<- updateThetaS(theta, S, sigma2,lam)
    beta<- invSoftThresholding(beta,theta,mu)
    
    sigma2<- updateSigma2(beta, tau, theta)
    beta<- updateBetaNotS(beta,tau,sigma2, mu)
    
    res<- randomWalkBeta(beta,sigma2,w, eps_beta,subset_size = 20)
    beta<- res[[1]]
    accept_rate_beta<- res[[2]]
    
    theta <- softThresholding(beta,mu)
    S<- (abs(beta)>mu)
    
    w<- updateW(S)
    
    
    if( step<burn){
      eps_beta = eps_beta* exp(accept_rate_beta - 0.234^(20/p))
      # print( paste(c("MH acceptance rates: ", accept_rate_mu_tilde, accept_rate_beta), collapse=" "))
      
    }else{
      trace_theta<- rbind(trace_theta, theta)
      trace_sigma2<- c(trace_sigma2, sigma2)
    }
    
    setTxtProgressBar(pb, step)
  }
  
  return (list(theta=trace_theta, sigma2 = trace_sigma2))
}