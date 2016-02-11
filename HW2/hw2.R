library(coda)
th = as.mcmc(rgamma(1e6,10,1012))
th.hpd = HPDinterval(th, prob = 0.90)
th.hpd
th.ci = qgamma(c(0.05,0.95),10,1012)

curve(dgamma(x, 10,1012), from=0, to=0.2, xlab=expression(theta),
      ylab=expression(paste("p(",theta,"|y)")))
lines(x=c(th.hpd[1], th.hpd[1], th.hpd[2], th.hpd[2]),
      y = qgamma(c(0, th.hpd[1], th.hpd[2], 0), 10, 1012), col="blue", lwd=2, lty=1)
lines(x=c(th.ci[1], th.ci[1], th.ci[2], th.ci[2]),
      y = qgamma(c(0, th.ci[1], th.ci[2], 0), 10, 1012), col="green", lwd=2, lty=1)
legend("topright", c("90% HPD interval", "90% credible interval"),
       col=c("blue", "green"), lty=c(1,2), lwd=c(2,2), bty="n")

