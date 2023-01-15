#############################################################
### Target trial emulation and causal inference practical ###
#############################################################

#################
### Section 2 ###
#################

set.seed(1)
lalonde <- read.table(file = "~/Downloads/IIAS-main/data/lalonde_nsw.csv") # change filepath accordingly

y <- lalonde$re78
y <- y/1000 # rescale y to $1000 units
t <- lalonde$treat
n <- length(y)

N_t <- sum(t)
N_c <- sum(1-t)

#############################
### Predictive resampling ###
#############################

M <- 3000 # number of forward simulation steps

tpred <- numeric(M)
ypred <- numeric(M)
N_t_cur <- N_t
N_c_cur <- N_c
sum_y1 <- sum(t*y)
sum_y0 <- sum((1-t)*y)

tpred <- numeric(M)
ypred <- numeric(M)

the_pred <- numeric(M-20)

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
  if (m > 20) {
    the_pred[m-20] <- sum(tpred[1:m]*ypred[1:m])/sum(tpred[1:m]) - sum((1-tpred[1:m])*ypred[1:m])/sum(1-tpred[1:m])
  }
  
}
plot(21:M, the_pred, type = "l", ylab = "ATE", xlab = "Forward step", ylim = c(-2,6))




