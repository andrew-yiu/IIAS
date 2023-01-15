#############################################################
### Target trial emulation and causal inference practical ###
#############################################################

###########################
### Section 1 solutions ###
###########################

set.seed(1)
lalonde <- read.table(file = "~/Downloads/lalonde_nsw.csv") # change filepath accordingly

## Question 1(i): compute the difference-in-means estimator

y <- lalonde$re78
t <- lalonde$treat
est <- sum(y*t)/sum(t) - sum(y*(1-t))/sum(1-t) # difference-in-means estimator
print(est)

## Question 1(iii): draw from the test statistic sampling distribution

rep <- 10000 # number of Monte Carlo samples to be drawn
n <- length(y)
est_vec <- numeric(rep) # create empty vector to store test statistic samples
for (r in 1:rep) {
  t_vec <- rbinom(n,1,prob=0.5) # generate random vector of treatment assignments
  y1_bar <- sum(y*t_vec)/sum(t_vec)
  y0_bar <- sum(y*(1-t_vec))/sum(1-t_vec)
  est_vec[r] <- y1_bar - y0_bar
}
hist(est_vec, xlab = "Test statistic", main = "Randomization distribution")
abline(v=est, col = "red")

## Question 1(iv): compute one-sided p-value estimate

length(which(est < est_vec))/rep # calculate proportion of samples less than the observed statistic

## Question 1(v)

beta <- 760
# fill in missing potential outcomes under the null
y1_vec <- y*t + (y + beta)*(1-t)
y0_vec <- y*(1-t) + (y - beta)*t

rep <- 100000 # increased no. of repetitions to add precision
n <- length(lalonde$treat)
est_vec_bet <- numeric(rep) 
for (r in 1:rep) {
  t_vec <- rbinom(n,1,prob=0.5)
  y1_bar <- sum(y1_vec*t_vec)/sum(t_vec)
  y0_bar <- sum(y0_vec*(1-t_vec))/sum(1-t_vec)
  est_vec_bet[r] <- y1_bar - y0_bar
}
# plot density estimates for sharp null and beta null
plot(density(est_vec), xlab = "Test statistic", main = "Randomization distribution")
lines(density(est_vec_bet), col = "green")

abline(v=est, col = "red")
bet0 <- expression(beta == 0)
bet760 <- expression(beta == 760)
legend(-2500, 6e-04, legend=c(bet0, bet760, "Observed"),
       col=c("black", "green", "red"), lty=1, cex=0.8)

# p-value estimate
length(which(est < est_vec_bet))/rep
