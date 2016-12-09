# Reference: Introducing Monte Carlo Methods With R - Springer

install.packages("mcsm")
library(mcsm)


# R basic functions -------------------------------------------------------

system.time(crossprod(1:10^6, 1:10^6)) # more efficient t(x)%*%y

# example of crossprod
matrix(c(1,2,3,4,5,6), nrow=2) %*% matrix(c(1,2,3,4,5,6), nrow=3)
crossprod(t(matrix(c(1,2,3,4,5,6), nrow=2)), matrix(c(1,2,3,4,5,6), nrow=3))
# end ex

# example of advanced matrix functions
chol(matrix) # upper triangle of Choleski decomposition
svd(matrix) # singular value decomposition
qr(matrix) # QR decomposition
ginv(matrix) # generalised inverse
backsolve(matrix) # inverse of upper diagonal triangular systems
forwardsolve(matrix) # inverse of lower diagonal triangular systems
chol2inv(matrix) # inverse of matrix m when provided by Choleski decomposition chol(m)
# end ex

array(1:50, c(3,5,3)) # create array

sample(1:12, 30, rep=T) # simulate 30 indep unif r.v. on {1,2,...,12}
sample(LETTERS[1:10], 30, rep=T) # simulate 30 indep unif r.v. on {a,b,...,j}

jitter(matrix(c(1,2,3,4,5,6), nrow=3)) # used to disrupt data with ties

# example of non-parametric kernel density estimates
Nit = c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,6,6,6)
AOB =c(4.26,4.15,4.68,6.08,5.87,6.92,6.87,6.25,
     + 6.84,6.34,6.56,6.52,7.39,7.38,7.74,7.76,8.14,7.22)
AOBm = tapply(AOB, Nit, mean)
Nitm = tapply(Nit, Nit, mean)
plot(Nit,AOB,xlim=c(0,6),ylim=c(min(AOB),max(AOB)),pch=19)

library(splines) #library for natural cubic splines

fitAOB=lm(AOBm ~ ns(Nitm, df=2)) #natural spline, fit depends on choice of ns
xmin=min(Nit); xmax=max(Nit)
lines(seq(xmin,xmax,.5), #fit to means
      predict(fitAOB,data.frame(Nitm=seq(xmin,xmax,.5))))

fitAOB2=loess(AOBm ~ Nitm,span = 1.25) #loess, fit depends on span
lines(seq(xmin,xmax,.5), #fit to means
        + predict(fitAOB2,data.frame(Nitm=seq(xmin,xmax,.5))))
#end ex

#bootstrap procedure use empirical distribution as a substitute for the true distribution
#to construct variance est. and CI

#bootstrap ex
x= seq(-3, 3, le=5) #equidispersed regressor
y = 2 + 4*x + rnorm(5) #simulated depedendent variable
fit = lm(y~x)
Rdata = fit$residuals

nBoot=2000 #set number of boostrap examples
B=array(0, dim=c(nBoot, 2)) # bootstrap array

for(i in 1:nBoot){ # bootstrap loop
  ystar = y + sample(Rdata, replace=T)
  Bfit = lm(ystar ~ x)
  B[i,] = Bfit$coefficients
}

colMeans(B) # show mean of the intercept and slope of the model above
sum((B[,1]-colMeans(B)[1])^2)/1999 #sample variance of the intercept
var(B[,1]) # same results as above
#end ex

#Cumulative sum for checking Monte Carlo convergence
x = rnorm(1)
for (t in 2:10^3)
  x = c(x, .09*x[t-1] + rnorm(1))

plot(x, type="l", xlab="time", ylab="x", lwd=2, lty=2,
     col="steelblue", ylim=range(cumsum(x))
)
lines(cumsum(x), lwd=2, col="orange3")
#end ex

#ex on polygon
par(mar=c(2,2,2,2))
x=matrix(0,ncol=100,nrow=10^4)
for (t in 2:10^4)
  x[t,]=x[t-1,]+rnorm(100)*10^(-2) # generate 100 Brownian motions in time [1, 10,000]

plot(seq(0,1,le=10^4),x[,1],ty="n",
     ylim=range(x),xlab="",ylab="")

polygon(c(1:10^4,10^4:1)/10^4, #x-coordinates
        c(apply(x,1,max),rev(apply(x,1,min))), #y-coordinates
        col="gold",bor=F) # range of 100 Brownian motions

polygon(c(1:10^4,10^4:1)/10^4,
        c(apply(x,1,quantile,.95),rev(apply(x,1,quantile,.05))), #90% confidence band
        col="brown",bor=F)
#end ex


# Random Variable Generation ----------------------------------------------

# check properties of Unif generator
Nsim=10^4 # number of random numbers
x=runif(Nsim)
x1=x[-Nsim]  # vectors to plot
x2=x[-1] # adjacent pairs
par(mfrow=c(1,3))
hist(x)
plot(x1,x2)
acf(x) # residual autocorrelation exisits in random numbers, common







