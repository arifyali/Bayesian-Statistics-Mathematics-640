##### Problem 1
library(MASS)
data(UScrime)
library(LearnBayes)

y = UScrime$y; n=length(y)
X = as.matrix(cbind(rep(1,n), UScrime[,-ncol(UScrime)]))

T = 1e4; k=dim(X)[2]; beta = matrix(NA, T, k)
g = n;  nu0 = 2;  s20 = 1

XtX.inv = solve(t(X)%*%X)
Sg = t(y)%*%(diag(1,n) - g/(g+1)*X%*%XtX.inv%*%t(X))%*%y
v = g/(g+1)*XtX.inv
m = v%*%t(X)%*%y

sigma2 = 1/rgamma(T, (nu0+n)/2, (nu0*s20+Sg)/2 )
for(t in 1:T) {
  beta[t,] = mvrnorm(1, m, sigma2[t]*v)
}

Crime.fit = lm(y~.,data= UScrime)
names(beta) = c("(Intercept)",names(UScrime)[-ncol(UScrime)])
Beta_means = apply(beta, 2, mean)
names(Beta_means) = c("(Intercept)",names(UScrime)[-ncol(UScrime)])
Beta_means
Crime.fit$coefficients

beta.ci = apply(beta, 2, quantile, c(0.025, 0.975))
colnames(beta.ci) = c("(Intercept)",names(UScrime)[-ncol(UScrime)])
t(beta.ci)
confint(Crime.fit)

Problem1b = function(seed = F, graph = 1){
  if(seed){
    set.seed(seed)
    }
  train = sample(nrow(UScrime), nrow(UScrime)*0.5)
  UScrime.train = UScrime[train,]
  UScrime.test = UScrime[-train,]
  Crimeb.fit = lm(y~.,data= UScrime.train)
  Crimeb.fit$coefficients
  Crimeb.fittedvalues = predict(Crimeb.fit, UScrime.test)
  mean((Crimeb.fittedvalues-UScrime.test$y)^2)

  y = UScrime.train$y; n=length(y)
  X = as.matrix(cbind(rep(1,n), UScrime.train[,-ncol(UScrime.train)]))

  T = 1e4; k=dim(X)[2]; beta = matrix(NA, T, k)
  g = n;  nu0 = 2;  s20 = 1

  XtX.inv = solve(t(X)%*%X)
  Sg = t(y)%*%(diag(1,n) - g/(g+1)*X%*%XtX.inv%*%t(X))%*%y
  v = g/(g+1)*XtX.inv
  m = v%*%t(X)%*%y

  sigma2 = 1/rgamma(T, (nu0+n)/2, (nu0*s20+Sg)/2 )
  for(t in 1:T) {
    beta[t,] = mvrnorm(1, m, sigma2[t]*v)
  }

  Beta_means = apply(beta, 2, mean)
  names(Beta_means) = c("(Intercept)",names(UScrime)[-ncol(UScrime)])
  Beta_means
  bayesian_predictions = (as.matrix(UScrime.test[,-ncol(UScrime.test)])%*%(Beta_means[-1])+Beta_means[1])[,1]
  
  mean((bayesian_predictions - UScrime.test$y)^2)
  
  if(graph != 0){
    par(mfrow=c(1,2))
    plot(UScrime.test$y, Crimeb.fittedvalues, main = "least-squares regression coefficients", xlab = "observations", ylab = expression(hat(y[i])))
    plot(UScrime.test$y, bayesian_predictions, main = "bayesian posterior means regression coefficients", xlab = "observations", ylab = expression(hat(y[i])))
    
  }
  
  results = list(Crimeb.fit$coefficients, mean((Crimeb.fittedvalues-UScrime.test$y)^2),Beta_means,mean((bayesian_predictions - UScrime.test$y)^2))
  names(results) = c("least-squares regression coefficients", "MSE", " posterior mean", "MSE  posterior mean")
  return(results)
}

Problem1b(1e8, 0)

Problem1c = matrix(data = NA, nrow = 100, ncol = 2)
names(Problem1c) = c("least-squares regression MSE", "MSE  posterior mean")
for(i in 1:100){
  Problem1ci = Problem1b()
  Problem1c[i,]= c(Problem1ci[[2]], Problem1ci[[4]])
}
apply(Problem1c,2,mean, na.rm = T)

#### Problem 2
msparrownest <- read.table("~/Documents/Georgetown/Bayesian-Statistics-Mathematics-640/HW5/msparrownest.dat", quote="\"", comment.char="")
