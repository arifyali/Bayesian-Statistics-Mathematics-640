# Homework 5
# Math 640, Georgetown
# Author: Alex Perrone

library("MCMCpack")
library("MASS")
library("coda")
library("reshape")
library("ggplot2")
library("data.table")
library("xtable")
library("magrittr")

options(digits=10)

################## 1. ##################
data("UScrime")  # requires MASS package

# 1(a)
bayes_fit <- function(dat){
  y <- dat$y
  n <- length(y)
  intercept <- rep(1, n)
  X <- as.matrix(cbind(intercept, dat[ , 1:15]))
  num_predictors <- ncol(X)
  
  # Hyperparameters for g prior.
  g <- n
  beta0 <- as.matrix(rep(0, num_predictors))  # zero column-vector
  sigma0_squared <- 1
  v0 <- 2
  
  # Compute quantities for posterior distribution of sigma^2.
  XtX_inv <- solve(t(X)%*%X)
  S_g_squared <- t(y) %*% 
    (diag(n) - (g / (g + 1)) * X %*% XtX_inv %*% t(X)) %*% 
    y
  
  # Draw posterior samples for sigma^2 | y, X. 
  niter <- 1e4
  phi <- rgamma(n = niter, 
                shape = 0.5 * (v0 + n), 
                rate = 0.5 * (v0 * sigma0_squared + S_g_squared))
  sigma2 <- 1 / phi
  
  # Compute quantities for posterior distribution of beta | sigma^2, y, X. 
  m <- g / (g + 1) * (beta0 / g + XtX_inv %*% t(X) %*% y)
  beta <- matrix(NA, niter, num_predictors)
  for (i in 1:niter){
    V <- g / (g + 1) * sigma2[i] * XtX_inv
    beta[i, ] <- mvrnorm(1, m, V)
  }
  beta <- data.table(beta)
  setnames(beta, attr(X, "dimnames")[[2]])  # extract col names from design matrix
  return(list(sigma2=sigma2, beta=beta))
}

# Trace plot of each variable. 
plot_bayes <- function(beta, plot=FALSE){
  melt_beta <- melt(beta)
  melt_beta[ , index := 1:.N, by=variable]
  if (plot){
    p <- ggplot(melt_beta, aes(x=index, y=value)) + 
      geom_line(size=0.05) + 
      facet_wrap(~ variable, scales = "free_y") + 
      labs(title="Trace Plot for Beta") + 
      xlab("Iteration") + 
      ylab("") + 
      theme(plot.margin = unit(c(0.5, 1, .5, .5), "cm"))
    ggsave("1a-trace.jpeg", plot = p, width=8, height=6)
    
    # Histogram of estimates. 
    p <- ggplot(melt_beta, aes(value)) + 
      geom_histogram(aes(y = ..density..)) + 
      facet_wrap(~ variable, scales = "free") + 
      labs(title="Histogram of Posterior Estimates for Beta") + 
      xlab("") + 
      ylab("Density") + 
      theme(plot.margin = unit(c(0.5, 1, .5, .5), "cm"))
    ggsave("1a-hist.png", plot = p, width=8, height=6)
  }
  
  # Posterior means. 
  mean_and_quantile <- function(x){
    list(mean=mean(x), 
         `95%-credible interval`=paste0("(", 
                                        round(quantile(x, 0.025), 3), 
                                        ", ", 
                                        round(quantile(x, 0.975), 3), 
                                        ")"))
  }
  posterior_beta <- melt_beta[ , mean_and_quantile(value), by=variable]
  latex <- data.frame(posterior_beta) %>% xtable(digits=3)
  print(latex, include.rownames=FALSE)
}

beta <- bayes_fit(UScrime)$beta
plot_bayes(beta, plot=TRUE)

# Frequentist least-squares estimate. 
lm1 <- lm(y ~ X)
summary(lm1) %>% xtable(digits=3) %>% print()

# 1(b)
# i. Least-squares estimates. 
train_test_split <- function(dat){
  n <- nrow(dat)
  index <- sample(n)  # permute the row numbers
  train_index <- index[1:floor(n / 2)]  # first half
  test_index <- index[(floor(n / 2) + 1):n]  # second half
  train <- dat[train_index, ]
  test <- dat[test_index, ]
  test <- test[order(test$y), ]
  return(list(train=train, test=test))
}

fit_lm <- function(train, test, plot=FALSE, plot_name=""){
  # Fit the lm. 
  lm_mod <- lm(y ~ ., data=train)
  
  # Generate predictions. 
  test$predicted <- predict.lm(lm_mod, newdata=test)
  mse <- mean((test$predicted - test$y)^2)
  
  # Plot, if desired. 
  if (plot){
    lm_mod %>% xtable(digits=3) %>% print()
    print(paste("MSE = ", mse))
    print(paste("RMSE = ", sqrt(mse)))
    plot_predicted_observed(test, plot_name)
  }
  return(mse)
}

plot_predicted_observed <- function(test, plot_name, credible=FALSE){
  # Plot predicted and observed. 
  test$index <- 1:nrow(test)
  test_plot_data <- melt(test[ , c("y", "predicted", "index")], 
                         id.vars="index")
  p <- ggplot(test_plot_data, aes(index, value, group=variable, 
                                  colour=variable, shape=variable)) + 
    geom_point(size=3) +
    labs(title="Predicted and Observed y") + 
    xlab("Index") + 
    ylab("y") + 
    theme(legend.title=element_blank())
  if (credible){
    mm <- melt(test[ , c("y", "predicted", "credible_interval_low", 
                         "credible_interval_high", "index")],
               id.vars=c("index", "credible_interval_low", 
                         "credible_interval_high"))
    p <- p + geom_errorbar(data=mm, aes(ymin=credible_interval_low, 
                                        ymax=credible_interval_high))
  } else {
    p <- p + geom_line() 
  }
  ggsave(plot_name, plot = p, width=8, height=6)
}

set.seed(1)  # make train/test split reproducible
split <- train_test_split(UScrime)
train <- split$train
test <- split$test
fit_lm(train, test, plot=TRUE, plot_name="1bi.png")

# Bayesian estimates. 
predict_bayes <- function(beta_train, test, plot=FALSE, plot_name="", 
                          credible=FALSE){
  beta_train <- as.matrix(beta_train)
  intercept <- rep(1, nrow(test))
  X <- as.matrix(cbind(intercept, test[ , 1:15]))
  res <- X %*% t(beta_train)
  test$predicted <- rowMeans(res)
  test$credible_interval_low <- apply(res, 1, quantile, 0.025)
  test$credible_interval_high <- apply(res, 1, quantile, 0.975)
  mse <- mean((test$predicted - test$y)^2)
  
  # Plot. 
  if (plot){
    print(mse)
    print(paste("MSE = ", mse))
    print(paste("RMSE = ", sqrt(mse)))
    plot_predicted_observed(test, plot_name, credible)  # turn credible = TRUE to show intervals
  }
  return(mse)
}
beta_train <- bayes_fit(train)$beta
plot_bayes(data.table(beta_train), plot=FALSE)  # do not make plots, but show table
mse <- predict_bayes(beta_train, test, plot=TRUE, plot_name="1bii.png")
mse <- predict_bayes(beta_train, test, plot=TRUE, 
                     plot_name="1bii-with-intervals.png", credible=TRUE)

# 1(c)
niter <- 100
mse <- data.frame(least_squares_mse=rep(NA, niter), 
                  bayesian_mse=rep(NA, niter))
for (i in 1:niter){
  print(i)
  split <- train_test_split(UScrime)
  train <- split$train
  test <- split$test
  mse$least_squares_mse[i] <- fit_lm(train, test, plot=FALSE)
  beta_train <- bayes_fit(train)$beta
  mse$bayesian_mse[i] <- predict_bayes(beta_train, test, plot=FALSE)
}
apply(mse, 2, mean)


################## 2. ##################
sparrow <- read.table("msparrownest.dat")
names(sparrow) <- c("nest", "wingspan")
y <- sparrow$nest
n <- length(y)
X <- as.matrix(cbind(rep(1, n), sparrow$wingspan))
sigma2_hat <- n * mean(y) * (1 - mean(y))
XtX_inv <- solve(t(X)%*%X)
sigma2_proposal <- sigma2_hat * XtX_inv

# Log-likelihood for theta. 
log_likelihood <- function(theta, y, X){
  yhat <- X %*% theta
  sum(yhat * y - log(1 + exp(yhat)))
}

# The priors for theta on the log-scale. 
log_prior <- function(theta){
  log(dnorm(theta[1], mean = 0, sd = 10)) + 
    log(dnorm(theta[2], mean = 0, sd = 10))
}

# The log-probability of theta. 
log_prob <- function(theta, y, X){
  log_likelihood(theta, y, X) + log_prior(theta)
}

# Metropolis algorithm. 
niter <- 1e5
theta <- matrix(NA, niter, 2)
set.seed(1)
theta[1, ] <- c(0, 0)  # uninformative prior values
for (i in 2:niter){
  theta_proposed <- mvrnorm(n = 1, mu = theta[i - 1, ], 
                            Sigma = sigma2_proposal)
  log_r <- log_prob(theta_proposed, y, X) - log_prob(theta[i - 1, ], y, X)
  r <- exp(log_r)
  u <- runif(1)
  if (u < r){
    theta[i, ] <- theta_proposed
  } else {
    theta[i, ] <- theta[i - 1, ]
  }
}

# Posterior means (no burn-in). 
apply(theta, 2, mean)

# Discard burn-in. 
burnin <- niter * 0.10
theta <- theta[-(1:burnin), ]

jpeg("2b-trace.jpeg", height=400, width=800)
par(mfrow=c(1,2))
plot(theta[, 1], type="l", main="Trace plot for alpha", 
     ylab="alpha", xlab="Iteration")
plot(theta[, 2], type="l", main="Trace plot for beta", 
     ylab="beta", xlab="Iteration")
dev.off()

# Effective sample size after burn-in removed. 
effectiveSize(as.mcmc(theta))

# Posterior means after burn-in removed. 
apply(theta, 2, mean)

# Histograms and posterior density estimate for alpha. 
jpeg("2c-alpha.jpeg", height=400, width=500)
hist(theta[, 1], breaks=30, xlim=c(-30, 30), probability=TRUE, 
     main="Prior and posterior for alpha", 
     xlab="alpha")
lines(density(theta[, 1]), col="orange", lwd=3)
# Prior density. 
x <- seq(-30, 30, length.out=1e4)
lines(x, dnorm(x, mean = 0, sd = 10), col="blue", lwd=3)
# Add legend. 
legend(x=6, y=0.09, c("Prior Density: N(0,100)", "Posterior Density"), 
       col=c("blue", "orange"), lwd=c(3, 3), bty="n")
dev.off()

# Histograms and posterior density estimate for beta.
jpeg("2c-beta.jpeg", height=400, width=500)
hist(theta[, 2], breaks=30, xlim=c(-30, 30), probability=TRUE, 
     main="Prior and posterior for beta", 
     xlab="beta")
lines(density(theta[, 2]), col="orange", lwd=3)
# Prior density. 
x <- seq(-10, 10, length.out=1e4)
curve(dnorm(x, mean = 0, sd = 10), col="blue", add=T)
# Add legend. 
legend(x=1.3, y=1, c("Prior Density: N(0,100)", "Posterior Density"), 
       col=c("blue", "orange"), lwd=c(3, 3), bty="n")
dev.off()

# 2(d)
# Setup. 
niter <- 1e3
start <- 10
end <- 16
wingspan_vec <- seq(start, end, length.out=niter)
predicted <- matrix(NA, niter, 4)
colnames(predicted) <- c("x", "f_ab_mean", 
                         "f_ab_credible_interval_low", 
                         "f_ab_credible_interval_high")
# Compute mean and credible interval. 
for (i in 1:length(wingspan_vec)){
  log_odds <- theta %*% c(1, wingspan_vec[i])  # alpha + beta * x
  f_ab <- exp(log_odds) / (1 + exp(log_odds))  # ratio, equal to probability of success
  predicted[i, ] <- c(wingspan_vec[i], mean(f_ab), 
                      quantile(f_ab, 0.025), quantile(f_ab, 0.975))
}
jpeg("2d-confidence-band.jpeg", width=400, height=320)
# Make a blank plot. 
plot(1, 1, type="n", xlim=c(start, end), ylim=c(0, 1), 
     main="Confidence band for f_ab", 
     xlab="Wingspan", ylab="f_ab( Wingspan )")
# https://stackoverflow.com/questions/14069629/plotting-confidence-intervals
polygon(c(rev(predicted[, 1]), predicted[, 1]), 
        c(rev(predicted[, 4]), predicted[, 3]), 
        col = "grey80", border = NA)
lines(predicted[, 1], predicted[, 2], lwd=3)  # mean of f_ab
lines(predicted[ , 1], predicted[, 3], lty="dashed", col="red")  # lower credible interval
lines(predicted[, 1], predicted[, 4], lty="dashed", col="red")   # upper credible interval
# Plot actual points from the data. 
points(sparrow$wingspan, sparrow$nest, pch=19, cex=0.8)
# Make a legend
legend(x=14.5, y=0.4, c("Observed", "Predicted"), bty="n",
       lty=c(0, 1), lwd=c(2, 2), col=c("black", "black"), 
       pch=c(19, NA), cex=0.8)
dev.off()
