#' @name l1ball
##' @title l1-ball prior for sparse regression
#' @description This package provides function for the l1-ball prior on high-dimensional regression. The main function, l1ball, yields posterior samples for linear regression.
#' @param y A data vector
#' @param X A design matrix
#' @param b_w The parameter in \eqn{Beta(1, p^{b_w})} for w
#' @param step Step size for the Markov Chain
#' @param burnin Number of burn-ins for the Markov Chain

#' @return trace returns a list of two component: trace_slab and trace_theta, containing all the posterior samples after burn-in.
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

#' @export l1ball



library(extraDistr)
Sample_a <- function(theta, lam, sigma2, slab_indicator, p){
	a = numeric(p)
	zero_idx = slab_indicator==0
	non_zero_idx = slab_indicator==1
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
SampleSlabIndicator <- function(p, mu, lam, a, theta, sigma2, slab_indicator, eps= .001){
	p_slab = -mu/lam/sqrt(sigma2)
	p_spike = log(1-exp(p_slab))

	sigma = sqrt(sigma2)
	w1 = p_spike + a_density(eps,theta,lam,sigma2)
	w2 = p_slab + a_density(a, theta, lam, sigma2)
	new_p = t(cbind(w1, w2)) + matrix(revd(p*2),nrow=2,ncol=p)
	slab =  apply(new_p, 2, which.max) ==2
	#print( sum(slab)/sum(slab_indicator))
	if (runif(1) < sum(slab)/sum(slab_indicator)){
		s <- slab
	}
	else{
		s <- slab_indicator
	}
	return(s)
}


SampleTheta <- function(X2, Xy, lam, a, sigma2, slab_indicator, p){
	a_star = numeric(p)
	idx = which(slab_indicator==1)
	a_star[idx] = a[idx]*lam[idx]^2
	theta = numeric(p)
	if (sum(slab_indicator)>0){
		a_starlam_sub_inv = diag(1/(a_star[idx]*lam[idx]^2))
		A = (X2[idx,idx]+a_starlam_sub_inv)
		theta0 = rnorm(sum(slab_indicator))
		LA = chol(A)
		#theta[idx] = solve(t(LA), theta0) * sqrt(sigma2) + solve(t(LA), solve(LA, Xy[idx]))
		theta[idx] = solve(t(LA), theta0) * sqrt(sigma2) + chol2inv(LA) %*% Xy[idx]

	}
	return(theta)
}

ExpCDF <- function(x,lam){
	return(1-exp(-x/lam))
}

ExpInvCDF <- function(y, lam){
	return(-log(1-y)*lam)
}

SampleT <- function(mu, lam, slab_indicator, theta, sigma2,p){
	m = lam*sqrt(sigma2)
	t = ExpInvCDF(ExpCDF(mu, m)*runif(p),m)-mu
	t[slab_indicator]=abs(theta)[slab_indicator]
	return(t)
}

SampleLam <- function(t, mu, sigma2, p){
	a=1
	b=.01
	beta = sum(t+mu)
	lam = 1/rgamma( n = p, shape = a+1, scale = 1/(beta/sqrt(sigma2))+b)
	return(lam)
}

SampleSigma2 <- function(X, y, theta, sigma2, t, mu, lam, n, p){
	beta = t+mu
	ig_m = lam* sqrt(sigma2)/beta
	a_beta = 1/rwald(p,ig_m,1)

	b1 = sum((y-X%*%theta)^2)/2 + sum(beta^2/lam^2/a_beta)/2
	a1 = n/2+ p/2 + .5 + .5
	return(1/rgamma(1, a1, 1/b1))
}

SampleMu <- function(p, lam, slab_indicator, sigma2, w, b_w, eps = .001){
	ComputeMu <- function(w, sigma2){
		a=1
		b=.01 * sqrt(sigma2)
		mu = (w ^ (-1.0/a)- 1)*b
		return(mu)
	}

	Compute_h_w <- function(w){
		mu = ComputeMu(w, sigma2)
		p_slab = -mu/lam/sqrt(sigma2)
		p_spike = log((1-exp(p_slab)))

		return(sum(p_slab*slab_indicator + p_spike*(1.0-slab_indicator))+(b_w-1)*log(1-w))
	}
	forward_lb = max(c(w-eps,0))
	forward_ub = min(c(w+eps, 1))
	forward_density = -log(forward_ub-forward_lb)
	w_new = runif(1,forward_lb, forward_ub)
	backward_lb = max(c(w_new - eps, 0))
	backward_ub = min(c(w_new+eps, 1))
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



l1ball <- function(y, X, b_w = 1.0, steps = 3000, burnin=1000){
	n = nrow(X)
	p = ncol(X)
	trace_theta = numeric()
	trace_slab = numeric()
	X2 = t(X) %*% X
	Xy = t(X) %*% y
	theta = solve(X2+diag(p), Xy)
	sigma2 = sum((y-X%*%theta)^2)/n
	mu = runif(1)
	lam = runif(p)*10
	slab_indicator = (abs(theta)>.05)
	w = 1./p

	trace_slab = matrix(0, nrow = step, ncol = p)
	trace_theta = matrix(0, nrow = step, ncol = p)
	for (k in 1:steps){
		if (k%%100==0){
			print(k)}
		a = Sample_a(theta, lam, sigma2, slab_indicator, p)
		slab_indicator = SampleSlabIndicator(p, mu, lam, a, theta, sigma2, slab_indicator, eps= .001)
		wmu = SampleMu(p, lam, slab_indicator, sigma2, w, b_w, eps = .001)
		w = wmu[1]
		mu = wmu[2]
		t = SampleT(mu, lam, slab_indicator, theta, sigma2,p)
		theta = SampleTheta(X2, Xy, lam, a, sigma2, slab_indicator, p)
		if (k>burnin){
			trace_theta[k,] <- theta
			trace_slab[k,] <- slab_indicator
		}
		lam = SampleLam(t, mu, sigma2, p)
		sigma2 = SampleSigma2(X, y, theta, sigma2, t, mu, lam, n, p)

	}
	return(list(trace_theta=trace_theta, trace_slab=trace_slab))
}
