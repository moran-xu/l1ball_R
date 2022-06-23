library('statmod')
library('truncnorm')


softThresholding<- function(beta,mu){
  sign(beta)*(abs(beta)-mu)*(abs(beta)>mu)
}


invSoftThresholding<- function(beta, theta, mu){
  S<- abs(theta)>0
  beta[S]<- sign(theta[S])*(abs(theta[S])+mu)
  beta
}

updateThetaS<- function(y,X,theta,tau, S, sigma2){
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

updateSigma2<- function(y,X,beta, tau, theta, n, p, a_sigma, b_sigma){
  alpha_sigma<- (n+p)/2+a_sigma
  beta_sigma<- (sum(beta*beta/tau)+ sum((y-X%*%theta)**2))/2+ b_sigma
  sigma2 = 1/rgamma(1,alpha_sigma,rate = beta_sigma)
  if(sigma2>10){
    sigma2=10 # truncate at 10 to prevent numeric issue
  }
  sigma2
}

updateLambda<-function(beta, sigma2, p, alpha, eta_tilde){
  rgamma( p, alpha+1, abs(beta)/sqrt(sigma2) + eta_tilde)
}

updateTau<- function(y,X,lam, sigma2, beta){
  p <- ncol(X)
  1/rinvgauss(p, mean= sqrt(lam**2* sigma2/(beta**2)), shape =lam**2)
}

updateBetaNotS<- function(beta, tau ,sigma2, S, mu){
  if(sum(!S)>0){
    beta[!S]<- rtruncnorm(sum(!S), a = -mu,b=mu, mean = 0, sd=sqrt(tau[!S]*sigma2))
  }
  beta
}

logPost<- function(beta, mu, sigma2, b_w, alpha, eta_tilde, y, X){
  theta<- softThresholding(beta, mu)
  part1<-  - sum((y-X%*%theta)**2)/2/sigma2
  part2<- - (alpha+1) * sum(log(1+ abs(beta)/eta_tilde))
  w<- (1+mu/eta_tilde)**(-alpha)
  part3 <- (b_w-1)*log(1-w)
  
  if(sum(abs(theta))==0){return(-Inf)}
  else{
    return(part1+ part2+ part3)
  }
}


randomWalkBeta<- function(y, X, beta, mu, b_w, sigma2, alpha, eta_tilde, eps=0.1,subset_size= 100){
  p <- ncol(X)
  accept =0
  mh_count =0
  
  cur_beta<- beta
  curLogPost<- logPost(beta, mu, sigma2, b_w, alpha, eta_tilde, y, X)
  
  full_idx<- sample(c(1:p))
  
  
  idx0<- 1
  idx1<- idx0 + subset_size-1
  
  while(idx0<p){
    if(idx1>p)idx1=p
    beta[full_idx[idx0:idx1] ]<- beta[full_idx[idx0:idx1]]+ runif(idx1-idx0+1,-1,1)*eps
    propLogPost<- logPost(beta, mu, sigma2, b_w,alpha, eta_tilde, y, X)
    
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


randomWalkMuTilde<- function(y, X, beta, mu_tilde, b_w, sigma2, alpha, eta_tilde, eps=0.1, tries= 20){
  p <- ncol(X)
  accept =0
  mh_count =0
  
  cur_mu_tilde <- mu_tilde
  curLogPost<- logPost(beta, abs(mu_tilde),sigma2, b_w, alpha, eta_tilde, y, X)
  
  i = 0
  while( i < tries){
    mu_tilde = mu_tilde + runif(1,-1,1)*eps
    propLogPost<- logPost(beta,abs(mu_tilde),sigma2, b_w, alpha, eta_tilde, y, X)
    mh_count = mh_count+1
    
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

l1ball <- function(y, X, b_w = 10, steps = 2000, burnin=500, alpha = .1, eta_tilde = 10., a_sigma = 2, b_sigma = 1, eps_beta = .1, eps_mu_tilde = .1, mu_tilde = .1, theta_initial = FALSE){
  n = nrow(X)
  p = ncol(X)
  trace_theta<- matrix(0, 0, p)
  trace_w <- numeric(steps)
  #### Initialization
  if (length(theta_initial)>1) { 
    theta0 <- theta_initial}
  else {theta0 <- lm.fit(X,y)$fitted.values}
  beta<- runif(p)
  tau = rep(1,p)
  mu<- 0.2
  mu_tilde <- 0.1
  sigma2 <- 0.1
  theta <- theta0
  S<- theta!=0
  beta <- invSoftThresholding(beta, theta, mu)
  beta[!S]= mu/2
  for(step in c(1:steps)){
    theta<- updateThetaS(y,X,theta,tau,S,sigma2)
    beta<- invSoftThresholding(beta,theta,mu)
    sigma2<- updateSigma2(y,X,beta, tau, theta, n, p, a_sigma, b_sigma)
    lam <- updateLambda(beta,sigma2, p, alpha, eta_tilde)
    tau<- updateTau(y,X,lam,sigma2,beta)
    beta<- updateBetaNotS(beta, tau, sigma2, S, mu)
    
    res<- randomWalkBeta(y, X, beta, mu, b_w, sigma2, alpha, eta_tilde, eps_beta, subset_size = 20)
    beta<- res[[1]]
    accept_rate_beta<- res[[2]]
    
    res<- randomWalkMuTilde(y, X, beta, mu_tilde, b_w, sigma2, alpha, eta_tilde, eps_mu_tilde, tries = 20)
    mu_tilde = res[[1]]
    accept_rate_mu_tilde = res[[2]]
    
    mu = abs(mu_tilde)
    
    theta <- softThresholding(beta,mu)
    S<- (abs(beta)>mu)
    trace_w[step] <- mu
    if(step < burnin){
      eps_beta = eps_beta* exp(accept_rate_beta - 0.234)
      eps_mu_tilde = eps_mu_tilde* exp(accept_rate_mu_tilde - 0.234)
    }else{
      trace_theta<- rbind(trace_theta, theta)
    }
    if (step %% 100 ==0)
      print( c('accept rate of mu_tilde:', accept_rate_mu_tilde, 'accept rate of beta:', accept_rate_beta))
  }
  list(trace_theta=trace_theta, trace_w = trace_w)
}

p <- 200
X <- diag(p)
w0 <- c(rnorm(5)/100+5,rep(0,p-5))
y <- X%*%w0+rnorm(p)
trace <- l1ball(y,X, b_w = 200., 2000, 1000, alpha = .1, eta_tilde =100.)
trace$trace_w

require("ggplot2")
require("reshape2")

df <- melt(trace$trace_theta[,])
colnames(df)<- c("name","index","value")
df$index<- as.factor(df$index)
# Basic boxplot

ggplot(df, aes(x=index, y=value)) +
  geom_boxplot(outlier.shape = NA) 
