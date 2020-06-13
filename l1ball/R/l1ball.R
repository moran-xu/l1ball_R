#' @name l1ball
##' @title l1-ball prior for sparse regression
#' @description This package provides function for the l1-ball prior on high-dimensional regression. The main function, l1ball, yields posterior samples for linear regression.
#' @param y A data vector
#' @param X A design matrix
#' @param b_w The parameter in \eqn{Beta(1, p^{b_w})} for w
#' @param step Step size for the Markov Chain
#' @param burnin Number of burn-ins for the Markov Chain

#' @return trace returns a list of two component: trace_NonZero and trace_theta, containing all the posterior samples after burn-in.
#' @import extraDistr
#' @importFrom extraDistr rwald
#' @import EnvStats
#' @importFrom EnvStats revd
#' @examples
#' p=100
#' X <- diag(p)
#' d =5
#' w0 <- c(rep(0, p-d), rnorm(d)+5)
#' y = X%*% w0
#' trace <- l1ball(y,X,1.5,3000,1000)
#'
#' @export l1ball



library(extraDistr)
Sample_a <- function(theta, lam, sigma2, NonZero_indicator, p){
	a = numeric(p)
	zero_idx = NonZero_indicator==0
	non_zero_idx = NonZero_indicator==1
	if (sum(non_zero_idx)>0){
		ig_m = lam[non_zero_idx] * sqrt(sigma2) / abs(theta[non_zero_idx])
		a[non_zero_idx] = 1/rwald(sum(non_zero_idx), ig_m, 1)
	}
	a[zero_idx] = rgamma(sum(zero_idx), shape = .5, scale = 2)
	return(a)
}

a_density <- function(a, theta, lam, sigma2){
	return(-log(a)/2 - theta ^2 /lam^2/2/sigma2/a-a/2)
}

library(EnvStats)

SampleNonZeroIndicator <- function(p, mu, lam, a, theta, sigma2, NonZero_indicator, eps= .001){
	p_NonZero = -mu/lam/sqrt(sigma2)
	p_Zero = log(1-exp(p_NonZero))

	sigma = sqrt(sigma2)
	w1 = p_Zero + a_density(eps,theta,lam,sigma2)
	w2 = p_NonZero + a_density(a, theta, lam, sigma2)
	new_p = t(cbind(w1, w2)) + matrix(revd(p*2),nrow=2,ncol=p)
	NonZero =  apply(new_p, 2, which.max) ==2
	#print( sum(NonZero)/sum(NonZero_indicator))
	if (runif(1) < sum(NonZero)/sum(NonZero_indicator)){
		s <- NonZero
	}
	else{
		s <- NonZero_indicator
	}
	return(s)
}


chol_ <- function(x){
  if(length(x)>0){
    return (chol(x))
    }
  else{
    return(sqrt(x))
  }

}


SampleTheta <- function(X2, Xy, lam, a, sigma2, NonZero_indicator, p){
	a_star = numeric(p)
	idx = which(NonZero_indicator==1)
	a_star[idx] = a[idx]*lam[idx]^2
	theta = numeric(p)
	if (sum(NonZero_indicator)>0){
		a_starlam_sub_inv = diag(1/(a_star[idx]*lam[idx]^2) + 1E-14) # add 1E-14 to make the cholesky work
		A = (X2[idx,idx]+a_starlam_sub_inv)
		theta0 = rnorm(sum(NonZero_indicator))
		LA = chol_(A)
		# theta[idx] = solve(t(LA), theta0) * sqrt(sigma2) + solve(t(LA), solve(LA, Xy[idx]))
		theta[idx] = solve(t(LA), theta0) * sqrt(sigma2) + chol2inv((LA)) %*% Xy[idx]

	}
	return(theta)
}

ExpCDF <- function(x,lam){
	return(1-exp(-x/lam))
}

ExpInvCDF <- function(y, lam){
	return(-log(1-y)*lam)
}

SampleT <- function(mu, lam, NonZero_indicator, theta, sigma2,p){
	m = lam*sqrt(sigma2)
	t = ExpInvCDF(ExpCDF(mu, m)*runif(p),m)-mu
	t[NonZero_indicator]=abs(theta)[NonZero_indicator]
	return(t)
}

SampleLam <- function(t, mu, sigma2, p){
	a=1
	b=1
	beta = sum(t+mu)
	lam = 1/rgamma( n = p, shape = a+1, scale = 1/(beta/sqrt(sigma2))+b)
	return(lam)
}

# SampleSigma2 <- function(X, y, theta, sigma2, t, mu, lam, n, p){
# 	beta = t+mu
# 	ig_m = lam* sqrt(sigma2)/beta
# 	a_beta = 1/rwald(p,ig_m,1)
# 	b1 = sum((y-X%*%theta)^2)/2 + sum(beta^2/lam^2/a_beta)/2
# 	a1 = n/2+ p/2 + .5 + .5
# 	return(1/rgamma(1, a1, 1/b1))
# }

SampleSigma2 <- function(X, y, theta, sigma2, t, mu, lam, n, p, eps_change = 1, ub= Inf){

  beta = t+mu

  ss2 = sum((y-X%*%theta)^2)
  s_beta = sum(abs(beta)/lam)

  Compute_l_sigma2<- function(sigma2){
    # using exponential prior on sigma \sigma \sim Exp(p) to prevent sigma2 increases too much
    - ((n+p)/2) * log(sigma2) - ss2/sigma2/2 - s_beta/sqrt(sigma2)    - p*sigma2
  }

  forward_lb = max(c(sigma2-eps_change,0))
  forward_ub = min( c(sigma2+eps_change, ub))

  forward_density = -log(forward_ub-forward_lb)
  sigma2_new = runif(1,forward_lb, forward_ub)

  backward_lb = max(c(sigma2_new - eps_change, 0))
  backward_ub = min(c(sigma2_new+eps_change,ub))
  backward_density = -log(backward_ub-backward_lb)

  if (log(runif(1))< Compute_l_sigma2(sigma2_new)+backward_density- Compute_l_sigma2(sigma2)-forward_density){
    sigma2 = sigma2_new
    accept =1
  }
  else{
    accept = 0
  }

  return( sigma2)

}

SampleMu <- function(p, lam, NonZero_indicator, sigma2, w, b_w, eps_change = 1E-2){

	ComputeMu <- function(w, sigma2){
		a=1
		b= sqrt(sigma2)
		mu = (w ^ (-1.0/a)- 1)*b
		return(mu)
	}

	Compute_h_w <- function(w){
		mu = ComputeMu(w, sigma2)
		p_NonZero = -mu/lam/sqrt(sigma2)
		p_Zero = log((1-exp(p_NonZero)))

		return(sum(p_NonZero*NonZero_indicator + p_Zero*(1.0-NonZero_indicator))+(p^b_w-1)*log(1-w))
	}
	forward_lb = max(c(w-eps_change,0))
	forward_ub = min(c(w+eps_change, 1))
	forward_density = -log(forward_ub-forward_lb)
	w_new = runif(1,forward_lb, forward_ub)
	backward_lb = max(c(w_new - eps_change, 0))
	backward_ub = min(c(w_new+eps_change, 1))
	backward_density = -log(backward_ub-backward_lb)

	if (log(runif(1))<Compute_h_w(w_new)+backward_density- Compute_h_w(w)-forward_density){
		w = w_new
		accept =1
	}
	else{
		accept = 0
	}
	mu = ComputeMu(w, sigma2)
	return(c(w, mu))
}



l1ball <- function(y, X, b_w = 1.0, steps = 3000, burnin=1000, eps=1E-8){
	n = nrow(X)
	p = ncol(X)


	trace_Sigma2 = numeric()
	trace_NonZero = matrix(0, nrow = steps, ncol = p)
	trace_theta = matrix(0, nrow = steps, ncol = p)
	trace_Lam = matrix(0, nrow = steps, ncol = p)

	X2 = t(X) %*% X
	Xy = t(X) %*% y
	theta = rnorm(p)#solve(X2+diag(p), Xy)
	sigma2 = 1 #sum((y-X%*%theta)^2)/n
	w = 1./p

	mu = quantile( abs(theta),1-w) #runif(1)
	lam = rep(1, p)
	NonZero_indicator = (abs(theta)>mu)


	for (k in 1:(steps+burnin)){
		if (k%%100==0){
			print(k)
		  # print(sigma2)
		  }
		a = Sample_a(theta, lam, sigma2, NonZero_indicator, p)
		NonZero_indicator = SampleNonZeroIndicator(p, mu, lam, a, theta, sigma2, NonZero_indicator, eps= eps)
		wmu = SampleMu(p, lam, NonZero_indicator, sigma2, w, b_w, eps_change = 1E-2)
		w = wmu[1]
		mu = wmu[2]
		t = SampleT(mu, lam, NonZero_indicator, theta, sigma2,p)
		theta = SampleTheta(X2, Xy, lam, a, sigma2, NonZero_indicator, p)

		lam = SampleLam(t, mu, sigma2, p)
		#sigma2 = 0.1**2#
		sigma2 = SampleSigma2(X, y, theta, sigma2, t, mu, lam, n, p)

		if (k>burnin){
		  idx = k-burnin
			trace_theta[idx,] <- theta
			trace_NonZero[idx,] <- NonZero_indicator
			trace_Sigma2[idx] <- sigma2
			trace_Lam[idx,] <- c(lam)

		}



	}
	return(list(trace_theta=trace_theta, trace_NonZero = trace_NonZero, trace_Sigma2= trace_Sigma2, trace_Lam= trace_Lam ))
}
