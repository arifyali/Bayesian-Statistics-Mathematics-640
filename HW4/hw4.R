U = runif(1e6,0,1)

X = log(U/(1-U))
hist(x)
curve(exp(x)/(1+exp(x))^2, add = T, col = "red")
mean(x<=-2)


f = function(x) {
  result = 30*x^2*(1-x)^2
    return(ifelse(result>=0 &result <=1,result, 0))
  }
g = function(x) return( dnorm(x,0.5,0.25) )
M = 1.2
x = seq(0, 1, length=1e6)
plot(x, M*g(x), type = "l", col = "blue", ylab="density")
lines(x, f(x), type = "l", col = "red")
legend("topright", legend = c("M*g(x)", "f(x)"), lty=c(1,1),
       col = c("blue", "red"), bty="n")
set.seed(1)
K=1e6
x.rs = accpt = rep(NA, K)
for(i in 1:K) {
  x = rnorm(1,0.5,0.25); u = runif(1)
  if(u < f(x)/(M*g(x))) {
    x.rs[i]=x 
    accpt[i]=1 }
  else accpt[i]=0
}
mean(x.rs, na.rm = T)

g = function(x) return( dbeta(x,2,2) )
M = 1.25
x = seq(0, 1, length=1e6)
plot(x, M*g(x), type = "l", col = "blue", ylab="density")
lines(x, f(x), type = "l", col = "red")
legend("topright", legend = c("M*g(x)", "f(x)"), lty=c(1,1),
       col = c("blue", "red"), bty="n")

K=1e6
x.rs = accpt = rep(NA, K)
for(i in 1:K) {
  x = rbeta(1,2,2); u = runif(1)
  if(u < f(x)/(M*g(x))) {
    x.rs[i]=x 
    accpt[i]=1 }
  else accpt[i]=0
}
mean(x.rs, na.rm = T)


samp.size=10000

x <- matrix(NA, samp.size, 2)
wts <- matrix(NA, samp.size, 2)

M = 1.2
x[,1] = rnorm(samp.size,0.5,0.25)
wts[,1] = f(x[,1])/(M*dnorm(x[,1],0.5,0.25))

M=1.25
x[,2] = rbeta(samp.size,2,2)
wts[,2] = f(x[,2])/(M*dbeta(x[,2],2,2))

mu.is = rep(NA, 0)
for(i in 1:2)
  mu.is[i] = sum((x[,i])*wts[,i])/sum(wts[,i])

mu.is




y <- read.table("~/Documents/Georgetown/Bayesian-Statistics-Mathematics-640/HW2/school.dat", quote="\"", comment.char="")[,1]
n = length(y)
T = 1e5
mu0 = 5
r20 = 0.5
alpha = 1
beta = 4

theta = matrix(NA, T, 2)
theta[1,] = c(30, 5)


for (i in 2:T) {
  theta[i,1] = rnorm(1, (sum(y)/theta[i-1,2]+mu0/r20)/(n/theta[i-1,2]+1/r20),
                     (n/theta[i-1,2]+1/r20)^-0.5)
  
  phi = rgamma(1, alpha+n/2,
               rate = 1/2*(beta+sum((y-theta[i,1])^2) ))
  theta[i,2] = 1/phi
}
nburn = 1e3
hist(theta[nburn:T,2], breaks=50,xlab=expression(sigma^2), prob=T, 
     col="orange", main=expression(paste("Posterior distribution of " ,sigma^2)))
lines(density(theta[,2]), lwd=2, col="black")



metro = function(sigma,theta0 = 0, T=1e5){
  theta = c(theta0)
  accept = c(0)
  print(sigma)
  accept = 0
  for(t in 2:T){
  theta_star = rnorm(1,theta[t-1], sigma)
  r = exp(-(theta_star^2-theta[t-1]^2)/2)
  u = runif(1)
  if(u<r){
    theta[t] = theta_star
    accept[t] = 1
  }else{
    theta[t] = theta[t-1]
    accept[t] = 0
  }
  }
  return(data.frame(theta,accept))
}

s1 = metro(0.5)
s2 = metro(2)

  par(mfrow=c(2,2))
  plot(s1[,1], type="l", xlab="iteration", ylab=expression(theta), 
        main=expression(paste("Trace plot of ",theta," for " ,sigma=0.5)))
  hist(s1[,1], breaks=50,xlab=expression(sigma), prob=T, 
        main=expression(paste("Plot of the Kernel Density of ",theta, " for ", sigma=0.5)))
  lines(density(s1[nburn:T,1]), lwd=2, col="red")
  
  plot(s2[,1], type="l", xlab="iteration", ylab=expression(theta), 
       main=expression(paste("Trace plot of ",theta," for " ,sigma=2)))
  hist(s2[,1], breaks=50,xlab=expression(sigma), prob=T, 
       main=expression(paste("Plot of the Kernel Density of ",theta, " for ", sigma=2)))
  lines(density(s2[nburn:T,1]), lwd=2, col="red")

nburn = 1e3
mean(s1[-(1:nburn),2]==1)
mean(s2[-(1:nburn),2]==1)

library(coda)
s1.mcmc = as.mcmc(s1[-(1:nburn),1])
s2.mcmc = as.mcmc(s2[-(1:nburn),1])

autocorr(s1.mcmc)
autocorr(s2.mcmc)

effectiveSize(s1.mcmc)
effectiveSize(s2.mcmc)
