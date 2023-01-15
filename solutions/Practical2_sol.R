#############################################################
### Target trial emulation and causal inference practical ###
#############################################################

###########################
### Section 2 solutions ###
###########################

set.seed(1)
lalonde <- read.table(file = "~/Downloads/lalonde_nsw.csv") # change filepath accordingly

y <- lalonde$re78
y <- y/1000 # rescale y to $1000 units
t <- lalonde$treat
n <- length(y)

N_t <- sum(t)
N_c <- sum(1-t)

## MLE's
muhat_t <- sum(t*y)/N_t
muhat_c <- sum((1-t)*y)/N_c
the_hat <- muhat_t - muhat_c

## posterior parameters
p_mean_c <- muhat_c * N_c*10000/(N_c*10000+25)
p_mean_t <- muhat_t * N_t*10000/(N_t*10000+64)
p_var_c <- (N_c/25+(1/10000))^(-1)
p_var_t <- (N_t/64+(1/10000))^(-1)

the_mean <- p_mean_t - p_mean_c
the_var <- p_var_c + p_var_t

## compute 95% central credible interval for theta
low <- round(qnorm(0.025, mean = the_mean, sd = sqrt(the_var)),3)*1000
upp <- round(qnorm(0.975, mean = the_mean, sd = sqrt(the_var)),3)*1000
print(paste0("Central 95% credible interval for the ATE: ($", low, ", $", upp, ")"))

the_samp <- rnorm(10000, mean = the_mean, sd = sqrt(the_var))
plot(density(the_samp), main = "Comparison of posterior densities", xlab = "ATE")

#############################
### Predictive resampling ###
#############################

M <- 3000

ypred <- numeric(0)
tpred <- numeric(0)
N_t_cur <- N_t
N_c_cur <- N_c
sum_y1 <- sum(t*y)
sum_y0 <- sum((1-t)*y)

tpred <- numeric(M)
ypred <- numeric(M)

the_pred <- numeric(M-20)

B <- 500
M <- 3000
theta_samp <- numeric(B)

for (b in 1:B) {
  ypred <- numeric(0)
  tpred <- numeric(0)
  N_t_cur <- N_t
  N_c_cur <- N_c
  sum_y1 <- sum(t*y)
  sum_y0 <- sum((1-t)*y)
  
  tpred <- numeric(M)
  ypred <- numeric(M)
  
  for (m in 1:M) {
    
    p_mean_c <- sum_y0*10000/(N_c_cur*10000+25)
    p_mean_t <- sum_y1*10000/(N_t_cur*10000+64)
    p_var_c <- (N_c_cur/25+(1/10000))^(-1)
    p_var_t <- (N_t_cur/64+(1/10000))^(-1)
    
    mu_c_draw <- rnorm(1, mean = p_mean_c, sd = sqrt(p_var_c))
    mu_t_draw <- rnorm(1, mean = p_mean_t, sd = sqrt(p_var_t))
    
    tnew <- rbinom(1, 1, 0.5)
    tpred[m] <- tnew
    N_t_cur <- N_t_cur + tnew
    N_c_cur <- N_c_cur + (1-tnew)
    
    if (tnew == 1) {
      ynew <- rnorm(1, mean = mu_t_draw, sd = 5)
      sum_y1 <- sum_y1 + ynew
    } else {
      ynew <- rnorm(1, mean = mu_c_draw, sd = 8)
      sum_y0 <- sum_y0 + ynew
    }
    ypred[m] <- ynew
    
  }
  theta_samp[b] <- sum(tpred*ypred)/sum(tpred) - sum((1-tpred)*ypred)/sum(1-tpred)
  
}

lines(density(theta_samp), col = "red")

##############################
### Inference for the SATE ###
##############################

B <- 5000 # number of posterior samples for SATE
SATE <- numeric(B)

## MLE's
muhat_t <- sum(t*y)/N_t
muhat_c <- sum((1-t)*y)/N_c
the_hat <- muhat_t - muhat_c

## posterior parameters
p_mean_c <- muhat_c * N_c*10000/(N_c*10000+25)
p_mean_t <- muhat_t * N_t*10000/(N_t*10000+64)
p_var_c <- (N_c/25+(1/10000))^(-1)
p_var_t <- (N_t/64+(1/10000))^(-1)

## draw posterior samples of SATE
for (b in 1:B) {
  # create data vectors for each potential outcome
  y1_vec <- t*y
  y0_vec <- (1-t)*y
  # draw mu parameters from their posteriors
  mu_c_draw <- rnorm(1, mean = p_mean_c, sd = sqrt(p_var_c))
  mu_t_draw <- rnorm(1, mean = p_mean_t, sd = sqrt(p_var_t))
  for (k in 1:n) {
    if (t[k] == 1) {
      # t=1 implies that the potential outcome for t=0 is missing
      y0_vec[k] <- rnorm(1, mean = mu_c_draw, sd = 5) 
    } else {
      # t=0 implies that the potential outcome for t=1 is missing
      y1_vec[k] <- rnorm(1, mean = mu_t_draw, sd = 8)
    }
  }
  SATE[b] <- mean(y1_vec) - mean(y0_vec) # compute SATE
  
}

## Compare posterior densities
plot(density(SATE), col = "green", xlim = c(-1,4), main = "Comparison of SATE and ATE posterior densities", xlab = "Treatment effect")
lines(density(the_samp))
lines(density(theta_samp), col = "red")
legend(-1, 1.2, legend=c("Exact ATE", "Pred. ATE", "SATE"),
       col=c("black", "red", "green"), lty=1, cex=0.8)
