### Exercise 1
# Part B
mu.c = rt(1e6, df=31)*sqrt(0.24^2/32)+1.013
mu.t = rt(1e6, df=35)*sqrt(0.2^2/36)+1.173
quantile(mu.t-mu.c, c(0.025,0.975))
hist(mu.t-mu.c, breaks = 50)
### Exercise 2
library(coda)
th = as.mcmc(rgamma(1e6,10,1012))
th.hpd = HPDinterval(th, prob = 0.90)
th.hpd
th.ci = qgamma(c(0.05,0.95),10,1012)

curve(dgamma(x, 10,1012), from=0, to=0.2, xlab=expression(theta),
      ylab=expression(paste("p(",theta,"|y)")))
lines(x=c(th.hpd[1], th.hpd[1], th.hpd[2], th.hpd[2]),
      y = dgamma(c(0, th.hpd[1], th.hpd[2], 0), 10, 1012), col="blue", lwd=2, lty=1)
lines(x=c(th.ci[1], th.ci[1], th.ci[2], th.ci[2]),
      y = dgamma(c(0, th.ci[1], th.ci[2], 0), 10, 1012), col="green", lwd=2, lty=1)
legend("topright", c("90% HPD interval", "90% credible interval"),
       col=c("blue", "green"), lty=c(1,2), lwd=c(2,2), bty="n")
z = seq(0,1000, 0.1)
density=11*(1012+95)^11/(1012+95+z)^12
plot(z,density, type = "l")
### Exercise 3
# Part A
school = read.table("~/Documents/Georgetown/Bayesian-Statistics-Mathematics-640/HW2/school.dat")[,1]
ybar = mean(school); s2 = var(school); n = length(school)
m0 = 0; k0 = 0.1; v0 = 10; s0 = 4
kn = k0 + n
mn = (k0*m0 + n*ybar)/kn
vn = v0 + n
sn = (v0*s0 + (n-1)*s2 + (k0*n/kn)*(ybar - m0)^2)/vn
phi  = rgamma(100000, vn/2, rate=vn/2*sn)
sigma2 = 1/phi
mu   = rnorm(100000, mn, sqrt(sigma2/kn))
# Part B
quantile(mu, c(0.025,0.975))
quantile(sigma2, c(0.025,0.975))
# Part C
mn + qt(c(0.025,0.975), df = vn) * sqrt(sn/kn)
# Part D
partd = rnorm(100000, mu, sqrt(sigma2))
# Part E
parte = mn + rt(100000, df = vn)*sqrt(sn*(1+1/kn))
qqplot(partd,parte)
