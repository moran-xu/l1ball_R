setwd("~/git/l1ball_R/")
rm(list=ls())


source("l1ball_anisotropic.R")
source("spike_gdp_slab.R")
source("aux.R")

n<- 200
p<- 200
d<- 5
sigma2_0<- 1
alpha_0<- 5 # true signal value

X<-  diag(n)
theta0<- rep(0,p)
theta0[1:d] <- rnorm(d, alpha_0, 0.01)
y<- rnorm(n, X%*%theta0, sqrt(sigma2_0))
# y[1]<- 5
y[1:d]


fit_l1ball<- l1ball(y,X, steps = 2000, burn=1000, 
                    alpha = 1, eta_tilde = 1,
                    a_w = 1, b_w = p) 
ts.plot(fit_l1ball$mu)

ts.plot(rowSums(fit_l1ball$theta!=0))

plotBoxplot(fit_l1ball,p_show = 200)

fit_sns<- spike_gdp_slab(y,X, steps = 5000, burn=1000,
                         alpha = 1, eta_tilde = 1,
                         a_w = 1, b_w = p
                         )
plotBoxplot(fit_sns,p_show = 100, FALSE)
ts.plot(rowSums(fit_sns$theta!=0))




