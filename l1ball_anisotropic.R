require(statmod)
require(truncnorm)

l1ball<- function(y,X, steps =1000,
                  alpha = 0.5, eta_tilde = 1,
                  a_sigma=2,b_sigma= 1, 
                  a_w = 1,b_w = p,
                  burn =500, theta_ini=NULL){
  
  softThresholding<- function(beta,mu){
    sign(beta)*(abs(beta)-mu)*(abs(beta)>mu)
  }
  
  
  invSoftThresholding<- function(beta,theta,mu){
    S<- abs(theta)>0
    beta[S]<- sign(theta[S])*(abs(theta[S])+mu)
    beta
  }
  
  updateThetaS<- function(theta, S, sigma2,lam){
    
    tau = rep(1,p)
    tau[S] = 1/rinvgauss(sum(S), mean= sqrt(lam[S]**2* sigma2/(theta[S]**2)), shape =lam[S]**2)
    
    
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
    alpha_sigma<- (n+p)/2+a_sigma
    beta_sigma<- (sum(beta*beta/tau)+ sum((y-X%*%theta)**2))/2+ b_sigma
    sigma2 = 1/rgamma(1,alpha_sigma,rate = beta_sigma)
    if(sigma2>10){
      sigma2=10 # truncate at 10 to prevent numeric issue
    }
    sigma2
  }
  
  updateLambda<-function(beta,sigma2){
    rgamma( p, alpha+1, abs(beta)/sqrt(sigma2) + eta_tilde)
  }
  
  updateTau<- function(lam,sigma2,beta){
    1/rinvgauss(p, mean= sqrt(lam**2* sigma2/(beta**2)), shape =lam**2)
  }
  
  updateBetaNotS<- function(beta,tau,sigma2, mu){
    if(sum(!S)>0){
      beta[!S]<- rtruncnorm(sum(!S), a = -mu,b=mu, mean = 0, sd=sqrt(tau[!S]*sigma2))
    }
    beta
  }
  
  getW<-function(mu, sigma2){
    (1+mu/eta_tilde/sqrt(sigma2))**(-alpha)
  }
  getMu<- function(w,sigma2){
    (w^(-1/alpha) -1) * eta_tilde * sqrt(sigma2)
  }
  
  
  logPost<- function(beta,mu,sigma2){
    theta<- softThresholding(beta,mu)
    part1<-  - sum((y-X%*%theta)**2)/2/sigma2
    part2<- - (alpha+1) * sum(log(1+ abs(beta)/eta_tilde/sqrt(sigma2)))
    w<- getW(mu,sigma2)
    part3 <- (a_w-1)*log(w)+(b_w-1)*log(1-w) 
    
    if(sum(abs(theta))==0){return(-Inf)}
    else{
      return(part1+ part2+ part3)
    }
  }
  
  
  randomWalkBeta<- function(beta,mu, eps=0.1,subset_size= 100){
    
    accept =0
    mh_count =0
    
    cur_beta<- beta
    curLogPost<- logPost(beta,mu,sigma2)
    
    full_idx<- sample(c(1:p))
    
    idx0<- 1
    idx1<- idx0 + subset_size-1
    
    while(idx0<p){
      if(idx1>p)idx1=p
      beta[full_idx[idx0:idx1] ]<- beta[full_idx[idx0:idx1]]+ rnorm(idx1-idx0+1,0,1)*eps  #+ runif(idx1-idx0+1,-1,1)*eps
      propLogPost<- logPost(beta,mu,sigma2)
      
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
  
  
  randomWalkMuTilde<- function(beta, mu_tilde, eps=0.1, tries= 20){
    
    accept =0
    mh_count =0
    
    curLogPost<- logPost(beta, abs(mu_tilde),sigma2)
    
    
    i = 0
    while( i < tries){
      
      S = abs(beta)> abs(mu_tilde)
      
      
      cur_mu_tilde <- mu_tilde
      
      curW<-  getW(abs(mu_tilde), sigma2 )
      
    
      newW<-  rbeta(1, curW/eps,   (1-curW)/eps) #rbeta(1, sum(S)+a_w, p-sum(S)+b_w)
      forward_prop = dbeta(newW, curW/eps, (1-curW)/eps, log = TRUE)
      mu_tilde <- getMu(newW, sigma2)
      
      newS = abs(beta)> abs(mu_tilde)
      backward_prop = dbeta(curW, newW/eps, (1-newW)/eps,log = TRUE)
      
      # mu_tilde = runif(1, max(abs(beta)[!S]), min(abs(beta)[S]))
      propLogPost<- logPost(beta,abs(mu_tilde),sigma2)
      mh_count = mh_count+1
      
      logRatio<-  propLogPost +  backward_prop - curLogPost - forward_prop 
      
      if (log(runif(1))< propLogPost -curLogPost){
        curLogPost = propLogPost
        cur_mu_tilde = mu_tilde
        accept = accept + 1
      }else{
        mu_tilde = cur_mu_tilde
      }
      
      i = i+1
    }
    return(list(mu_tilde, accept/mh_count))
  }
  
  
  
  # randomWalkMuTilde<- function(beta, mu_tilde, eps=0.1, tries= 20){
  #   
  #   accept =0
  #   mh_count =0
  #   
  #   curLogPost<- logPost(beta, abs(mu_tilde),sigma2)
  #   
  #   i = 0
  #   while( i < tries){
  #     
  #     S = abs(beta)>abs(mu_tilde)
  #     
  #     cur_mu_tilde <- mu_tilde
  #     
  #     # mu_tilde = runif(1, max(abs(beta)[!S]), min(abs(beta)[S]))
  #     mu_tilde = mu_tilde + runif(1,-1,1)*eps
  #     propLogPost<- logPost(beta,abs(mu_tilde),sigma2)
  #     mh_count = mh_count+1
  #     
  #     if (log(runif(1))< propLogPost -curLogPost){
  #       curLogPost = propLogPost
  #       cur_mu_tilde = mu_tilde
  #       accept = accept + 1
  #     }else{
  #       mu_tilde = cur_mu_tilde
  #     }
  #     
  #     i = i+1
  #   }
  #   return(list(mu_tilde, accept/mh_count))
  # }
  # 
  
  
  
  n<- nrow(X)
  p<- ncol(X)
  
  beta<- runif(p)
  
  
  sigma2<- 1
  
  # w= (1+mu/eta_tilde/sqrt(sigma2))**(-alpha)
  #mu<- (0.5^(-1/alpha) -1) * eta_tilde * sqrt(sigma2)
  
  mu<- 0.2
  
  mu_tilde = mu
  
  
  # theta<-  theta0 
  
  if(is.null(theta_ini)){
    theta<- solve(t(X)%*%X + diag(0.1,p), t(X)%*%y)
  }else{
    theta<- theta_ini  
  }
  S<- theta!=0
  beta<- invSoftThresholding(beta,theta,mu)
  beta[!S]= mu/2
  lam <- rep(1,p)
  
  tau<- rep(1,p)
  
  eps_beta= 0.1
  eps_mu_tilde = 0.1
  
  
  
  
  
  
  trace_theta<- matrix(0,0,p)
  trace_mu<- numeric(0)
  trace_sigma2<- numeric(0)
  
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = steps, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  subset_size = 20
  
  for(step in c(1:steps)){
    theta<- updateThetaS(theta,S,sigma2,lam)
    beta<- invSoftThresholding(beta,theta,mu)
    lam <- updateLambda(beta,sigma2)
    tau<- updateTau(lam,sigma2,beta)
    sigma2<- updateSigma2(beta, tau, theta)
    
    
    beta<- updateBetaNotS(beta, tau,sigma2, abs(mu_tilde))
    
    res<- randomWalkBeta(beta,mu,eps_beta,subset_size = subset_size)
    beta<- res[[1]]
    accept_rate_beta<- res[[2]]
    
    res<- randomWalkMuTilde(beta,mu_tilde, eps_mu_tilde, tries = 20)
    mu_tilde = res[[1]]
    accept_rate_mu_tilde = res[[2]]
    # accept_rate_mu_tilde=0.4
    
    mu = abs(mu_tilde)
    
    theta <- softThresholding(beta,mu)
    S<- (abs(beta)>mu)
    
    if( step<burn){
      eps_beta = eps_beta* exp(accept_rate_beta - 0.4^(subset_size/p))
      eps_mu_tilde = eps_mu_tilde* exp(accept_rate_mu_tilde - 0.4)
      
      #print( paste(c("MH acceptance rates: ", accept_rate_mu_tilde, accept_rate_beta), collapse=" "))
      
    }else{
      trace_theta<- rbind(trace_theta, theta)
      trace_mu<- c(trace_mu,mu)
      trace_sigma2<- c(trace_sigma2, sigma2)
    }
    
    setTxtProgressBar(pb, step)
  }
  
  return (list(theta=trace_theta, mu=trace_mu, sigma2 = trace_sigma2))
}