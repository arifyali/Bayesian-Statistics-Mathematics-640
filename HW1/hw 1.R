### 2c
n = 10
theta = c((1/(40^2)*180+n/(20^2)*150)/(1/(40^2)+n/(20^2)), 1/(1/40^2+n/20^2))
yhat = c((1/(40^2)*180+n/(20^2)*150)/(1/(40^2)+n/(20^2)), 1/(1/40^2+n/20^2)+20^2)
qnorm(c(0.025, 0.975), mean = theta[1], sd = sqrt(theta[2]))
qnorm(c(0.025, 0.975), mean = yhat[1], sd = sqrt(yhat[2]))

### 2d
n = 100
theta = c((1/(40^2)*180+n/(20^2)*150)/(1/(40^2)+n/(20^2)), 1/(1/40^2+n/20^2))
yhat = c((1/(40^2)*180+n/(20^2)*150)/(1/(40^2)+n/(20^2)), 1/(1/40^2+n/20^2)+20^2)
qnorm(c(0.025, 0.975), mean = theta[1], sd = sqrt(theta[2]))
qnorm(c(0.025, 0.975), mean = yhat[1], sd = sqrt(yhat[2]))

### 3c
r = 5
Y = c(7,10,5,8,6,12,6,9,7,5)
alpha = c(1,.5,5,1,5)
beta = c(1,.5,5,5,1)
par(mfrow = c(3,1))
for(i in 1:length(alpha)){
  theta = seq(0,1,0.01)
  credible_interval = qbeta(c(0.025,0.975), r*length(Y)+alpha[i], sum(Y)+beta[i])
  posterior_mean = (r*length(Y)+alpha[i])/(r*length(Y)+alpha[i]+sum(Y)+beta[i])
  plot(theta,  dbeta(theta, r*length(Y)+alpha[i], sum(Y)+beta[i]), 
       type = "l", col = "blue", xlab = expression(theta), 
       ylab = "Density",
       main = paste("alpha = ", alpha[i],"and", "beta = ", beta[i],
                    "- blue is posterior and red is prior"),
       sub = paste("The Credible Interval is", "(",credible_interval[1],",",
                   credible_interval[2], ")",
                   "and the posterior mean is", posterior_mean))
  lines(theta, dbeta(theta, alpha[i], beta[i]), col = "red")
}
