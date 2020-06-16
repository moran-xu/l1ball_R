rm(list=ls())

setwd("~/git/l1ball_R/")
source("l1ball/R/l1ball.R")

#Regression
n = 200
p = 500
X <- matrix(rnorm(n*p),n,p)
d =5
w0 <-   - c(rep(0, p-d), rnorm(d)*0.1 +2)
y = X%*% w0 + rnorm(n,0,.5)
trace <- l1ball(y, X, steps = 1000,burnin = 2000)

plot(trace$trace_Sigma2)

plot(w0)
plot(apply(trace$trace_theta,2,mean))

# plot(trace$trace_theta[,p-d+1])
# acf(trace$trace_theta[,p-d+1])



#change point detection

p = 500
n = 500
d = 5
idx = floor(seq(1,n,length.out=(d+2))[2:(d+1)])
w0 <- rep(0, p)
w0[idx] = -( rnorm(length(idx)) + 5)
X <- matrix(1, p, p)
X[upper.tri(X)] = 0
y = X%*%w0 + rnorm(p)*1
plot(y)
trace <- l1ball(y, X, steps = 1000,burnin = 2000)

theta_mean <- apply(trace$trace_theta,2,mean)

plot(y)
lines(X%*%theta_mean,col="red")

# plot(trace$trace_theta[,idx[1]])
# acf(trace$trace_theta[,idx[1]])



#Linear trend filter

n = 500
p = 500
d = 3

idx = floor(seq(1,n,length.out = (d+2))[2:(d+1)])

X= matrix(0,n,n)
for (i in c(1:n)){
  for (j in c(1:i)){
    X[i,j] = abs(j-i)+1
  }
}
X[,1]=1

# there are 4 line segments

idx= c(2,idx)

w0 = rep(0,p)
w0[idx] = (c(1,-2,1,-2))
psi = X%*%w0

y = psi + rnorm(n)*1
plot(y)

trace <- l1ball(y, X, steps = 2000,burnin = 2000)

theta_mean <- apply(trace$trace_theta,2,mean)


plot(theta_mean)
lines(w0)

plot(y)
lines(X%*%theta_mean,col="red")


