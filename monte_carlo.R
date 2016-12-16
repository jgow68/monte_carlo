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

# inverse transform method
Nsim=10^4 #number of random variables
U=runif(Nsim)
X=-log(U) #transforms of uniforms
Y=rexp(Nsim) #exponentials from R
par(mfrow=c(1,2)) #plots
hist(X,freq=F,main="Exp from Uniform")
hist(Y,freq=F,main="Exp from R")
rchisq(); rlogis(); rcauchy() # RNG for chi-squared, logistic, Cuchy distributions

# simulation of multivariate normal?
install.packages("mnormt")
library(mnormt)
rmnorm()
sadmvn(low=c(1,2,3),upp=c(10,11,12),mean=rep(0,3),var=B) # computation of probability of hypercubes?
#end ex

# another ex of multi-normal simulation (3x3 matrix)
Sigma=cov(matrix(rnorm(30),nrow=10)) 
A = t(chol(Sigma)) # Cholesky decomposition of Sigma to A%*%t(A)
x = A%*% rnorm(3) # 3-dim normal vector
apply(c, 1, function(x) A%*%rnorm(3)) # generate 100 vectors
#end ex

# discrete distributions

# binomial example ~ Bin(n,p)
# generate cdf for x=0,1,...,n
# generate U[0,1] then allocate random number k if P(x<k-1) < u < P(x<k)
# most observations will fall between interval mean +- 3 s.d.

# generate Poisson r.v. with value lambda
Nsim=10^4; lambda=100
spread=3*sqrt(lambda)
t=round(seq(max(0,lambda-spread),lambda+spread,1)) #integer values in the range around mean
prob=ppois(t, lambda)
X=rep(0,Nsim)
for (i in 1:Nsim){
  u=runif(1)
  X[i]=t[1]+sum(prob<u)} #doing so min(X)=min(t)
system.time(rpois(n=10^4, lambda=100)) #more efficient way of generating random numbers
#end ex

# Mixture representations
# example of negative binomial distribution ~ Neg(n,p)
Nsim=10^4; n=6; p=0.3;
y=rgamma(Nsim, n, rate=p/(1-p))
x=rpois(Nsim,y) # generate negative binomial r.v. from Poisson conditioned on Gamma distribution
hist(x, main="", freq=F, col="grey", breaks=40)
lines(1:50, dnbinom(1:50, n, p), lwd=2, col="sienna") # fit negative binomial to histogram
#end ex

# Accept-reject methods
# simulate beta~Be(alpha, beta)
optimise(f=function(x){dbeta(x, 2.7, 6.3)}, interval=c(0,1), maximum=T)$objective #maximum of the beta density in interval (0,1)
M=2.67 #upper bound
Nsim=2500; a=2.7; b=6.3
u=runif(Nsim, max=M) #uniform over (0, M) due to U x M
y=runif(Nsim) # generation from g
x=y[u<dbeta(y,a,b)] # u < f(y)/g(y) but g(y)=1 here

1/M #probability of acceptance of a given simulation
length(x)/Nsim # simulated no. of acceptances close to 1/M
M #is also the expcted no. of runs to get one acceptance, try to get small M
#end ex

# simulate beta r.v. with different method
x=NULL
while(length(x)<Nsim){
  y=runif(Nsim*M)
  x=c(x, y[runif(Nsim*M)*M < dbeta(y,a,b)])
}
x=x[1:Nsim]
mean(x)
#end ex

# difference between the 2 methods used above lies with runif
sim1 = runif(Nsim, max=M) #unif with interval (0,M)
sim2 = runif(Nsim*M) #unif(0,1) but with higher number of simulations
min(sim1); min(sim2)
max(sim1); max(sim2)
mean(sim1); mean(sim2)
#end ex

# simulate the beta function with Be(2,6)
optimise(f=function(x){dbeta(x,a,b)/dbeta(x,2,6)},
         maximum=T, interval=c(0,1))$objective #lower M value, hence higher acceptance rate
#end ex

# generate N(0,1) from a double-exponential with denstiy g(x|a) = (a/2)exp(-a*mod(x))
optimise(f=function(a){sqrt(2/pi)*(a^-1)*exp((a^2)/2)}, #optimise f(x)/g(x|a), f(x) is pdf of std normal
         interval=c(0,1), maximum=F) #minimum when a=1
#end ex

# ex on Weilbull(a,b) for Approx Bayesian Computation (ABC) Rejection Test
observedData = rweibull(20, 2, 5)
observedSummary = c(mean(observedData), sd(observedData))

model = function(par){
  simulatedData = rweibull(20, par[1,1], par[1,2])
  simulatedSummary = c(mean(simulatedData), sd(simulatedData))
  return(simulatedSummary)
}

n = 2*10^4
fit = data.frame(shape=runif(n, 0.01, 6), scale=runif(n, 0.01, 10),
                 summary1 = rep(NA, n), summary2 = rep(NA, n),
                 distance=rep(NA, n))

for (i in 1:n){ #for loop not the most efficient?
  prediction = model(fit[i, 1:2]) #simulate n Weibull with random parameters 
  deviation = sqrt(sum((prediction - observedSummary)^2)) #Euclidean distance, better if weigh outputs by variance
  fit[i, 3:5] = c(prediction, deviation)
}

# filter for simulated figures close to observed, distance < epsilon
plot(fit[fit[,5] < 1.5, 1:2], xlim=c(0.01, 6), ylim=c(0.01,10), col="lightgrey",
     main="Accepted parameters for different values of epsilon")
points(fit[fit[,5] < 1, 1:2], pch=18, col="grey")
points(fit[fit[,5] < 0.5, 1:2], pch=8, col="red")

legend("topright", c("<1.5", "<1", "<0.5"), 
       pch=c(1,18,8),
       col=c("lightgrey", "grey", "red")
       )

abline(v=2) #observed shape
abline(h=5) #observed scale
# accepted simulations should fall near the observed values

#end ex

# ex simulation on trigo function
x = seq(0,1, by=0.001)
y1 = (exp(-x^2/2)*(sin(6*x)^2 + 3*cos(x)^2*sin(4*x)^2 + 1)) #density that we want to simulate
y2 = exp(-x^2/2)/(sqrt(2*pi)) # we use simpler std normal density

optimise(f=function(x){(sin(6*x)^2 + 3*cos(x)^2*sin(4*x)^2 + 1)*(sqrt(2*pi))}, interval=c(0,1), maximum=T)
# bound M = 10.94

curve((exp(-x^2/2)*(sin(6*x)^2 + 3*cos(x)^2*sin(4*x)^2 + 1))) # ploting of original density

plot(x, y1, col="yellow", ylim=c(0,4))
points(x, y2*M, col="red") # std normal plot with bound M

y3 = y1/(y2*M) # accept-reject density
points(x, y3, col="grey")

Nsim=2500; M=10.94
sim1=NULL
while(length(sim1)<Nsim){
  y=rnorm(Nsim*M) # generate from g which is std Normal
  sim1=c(sim1, y[runif(Nsim*M)*M < (sqrt(2*pi)*(sin(6*y)^2 + 3*cos(y)^2*sin(4*y)^2 + 1))])
}
sim1=abs(sim1[1:Nsim]) # convert all values to positive

hist(sim1)
hist(y1)
par(mfrow=c(1,2))
# end ex