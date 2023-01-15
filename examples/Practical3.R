#############################################################
### Target trial emulation and causal inference practical ###
#############################################################

#################
### Section 3 ###
#################

nhanes <- read.table(file = "~/Downloads/nhanes_sub.csv")

y <- nhanes$years.lived.since.1971.up.to.1992
t <- nhanes$physically.active
x <- nhanes$age.at.interview

n <- length(y)

#########################
### Linear regression ###
#########################

B <- 200 ## number of posterior samples
M <- 2000 ## number of forward steps in predictive resampling

ate_samp <- numeric(B) ## create empty vector to store ATE samples

### compute initial parameter values
dat <- data.frame(x = x, y = y, t = t)
lin_fit <- lm(y ~ t + x, data = dat)
coef <- lin_fit$coefficients
w <- cbind(rep(1, length(x)), t, x)
s2 <- (1/(n-3))*sum((y-w %*% as.matrix(coef))^2)
V <- solve(t(w) %*% w)

for (b in 1:B) {
  N <- length(y) # N = current total size of data (observed + target trial)
  
  ### initialize parameter values from observed data
  ypred <- y
  s2pred <- s2
  Vpred <- V
  coefpred <- coef
  
  ### Bayesian bootstrap "shortcut"
  dir_wts <- rexp(n)
  dir_wts <- dir_wts/sum(dir_wts)
  x_count <- sample(1:n, size = M, replace = TRUE, prob = dir_wts)
  xpred <- x[x_count]
  # xpred <- rnorm(M, mean = 60, sd = 8) 
  xpred <- c(x, xpred)
  
  ### Generate all fully randomized treatment assignments
  tpred <- c(t, rbinom(M, 1, prob = 0.5))
  
  ### Combine t and x together into w
  wpred <- cbind(rep(1,n+M), tpred, xpred) # vector of 1's is for the intercept
  for (j in 1:M) {
    ### draw new y value from generalized t-distribution
    ynew <- rt(1,df = N-3)
    ynew <- sum(as.vector(wpred[n+j,])*coefpred) + ynew * sqrt(s2pred*(1+t(wpred[n+j,]) %*% Vpred %*% wpred[n+j,]))
    
    ### append new y value to the set of already observed outcomes
    ypred <- c(ypred, ynew)
    
    ### update parameters
    Vpred <- solve(crossprod(wpred[1:(n+j),]))
    N <- N+1
    coefpred <- tcrossprod(Vpred, wpred[1:(n+j),]) %*% as.matrix(ypred)
    s2pred <- (1/(N-3))*sum((ypred-wpred[1:(n+j),] %*% as.matrix(coefpred))^2)
    
  }
  ### compute ATE via difference-in-means estimator for the imputed target trial population
  ysim <- ypred[(n+1):(n+M)]
  tsim <- wpred[(n+1):(n+M), 2]
  ate_samp[b] <- mean(ysim[which(tsim==1)]) - mean(ysim[which(tsim==0)])
}
hist(ate_samp, main = "Posterior for the ATE", xlab = "ATE (years)", prob = TRUE)
lines(density(ate_samp), col = "red")
