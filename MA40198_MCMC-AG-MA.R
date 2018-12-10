rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")


rats$litter <- factor(rats$litter)

Z <- model.matrix(~ litter - 1, rats)
X <- model.matrix(~ 1 + rx, rats)
y <- rats$time



# Expressions for model including random effects
weib_re_prob <- deriv(
  expression(
    -log_sigma - eta + (1/exp(log_sigma) - 1)*(log(t) - eta) - (t/exp(eta))^(1/exp(log_sigma))
  ),
  namevec = "eta",
  function.arg = c("t", "eta", "log_sigma"),
  hessian = TRUE
)


weib_re_surv <- deriv(
  expression(
    -(t/exp(eta))^(1/exp(log_sigma))
  ),
  namevec = "eta",
  function.arg = c("t", "eta", "log_sigma"),
  hessian = TRUE
)




log_posterior <- function(theta, y, b, X, Z) {
  
  beta <- theta[1:2]
  log_sigma <- theta[3]
  sigma_b <- theta[4]
  
  eta <- X%*%beta + Z%*%b
  
  status <- X[, 2]
  
  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_re_prob(y, eta, log_sigma),
    (1-status) * weib_re_surv(y, eta, log_sigma)
  )
  
  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = sigma_b, log = TRUE))
  
  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb
  
  # Define log-prior (sum of log-priors for all parameters - here, just for log_sigma_b)
  log_prior <- dexp(sigma_b, rate = 5, log = TRUE)
  
  # Log-posterior is log-prior plus log-likelihood
  lf + log_prior
  
}


log_posterior(c(1, 1, 1, 1), rats$time, rep(0.1, ncol(Z)), X, Z)

set.seed(1);n.mc <- 100000;
theta <- c(1,1,1,1);b <- rep(0,50);
log_post_prop <- log_post <- log_posterior(theta, y, b, X, Z)
accept.b <- accept.th <- 0 ## acceptance counters
th <- matrix(0,4,n.mc); th[,1] <- theta; B <- matrix(0,50,n.mc) 
for (i in 2:n.mc) { ## Metropolis Hastings loop 
  ## update theta... 
  theta <- theta + rnorm(4)*c(0.01, 0.01,0.01,0.01)
  if (theta[4] > 0) log_post_prop <- log_posterior(theta, y, b, X, Z)
  if (theta[4] > 0 && runif(1) < exp(log_post_prop-log_post)) { 
    ## accept 
    accept.th <- accept.th + 1 
    log_post<-log_post_prop
  } else theta <- th[,i-1] ## reject 
  ## update random effects... 
  b <- b + rnorm(50)*0.01 
  log_post_prop <-log_posterior(theta, y, b, X, Z)
  if (runif(1) < exp(log_post_prop-log_post)) { ## accept 
    accept.b <- accept.b + 1 
    log_post<-log_post_prop
  } else b <- B[,i-1] ## reject
  th[,i] <- theta;B[,i] <- b 
} ## end of MH loop 
c(accept.b,accept.th)/n.mc ## acceptance rates

par(mfrow=c(2,2)) 
variable_names = c("beta0", "beta1", "log_sigma", "sigma_b")
for (j in 1:4) plot(th[j,],type="l", xlab=variable_names[j]) 




burn_in=10000
th_new=th[,burn_in:dim(th)[2]]
B_new = B[,burn_in:dim(B)[2]] #Get rid of 'burn-in' period

for (j in 1:4) plot(th_new[j,],type="l", xlab=variable_names[j]) 

par(mfrow=c(2,2)) 
for (j in 1:4) acf(th_new[j,]) #Slow mixing - very very high correlation


#Effective sample size, and checking for convergence using the Kolmogorov-Smirnov test. NOT WORKING - KS TEST SAYING NOT FROM SAME DISTRIBUTION.
n.eff <- c(0,0)
t.eff <- c(0,0)
par(mfrow=c(2,2))
for (j in 1:4) {
  acor <- acf(th_new[j,])
  t.eff[j] <- 2*sum(acor[[1]]) + 1
  n.eff[j] <- dim(th_new)[2]/t.eff[j]
}
n.eff<- as.integer(n.eff) #From Lec 9, this is the 'effective sample size' for the 4 parameters
t.eff<-as.integer(t.eff)


for (i in 1:4){
  indep_seq = seq(1, dim(th_new)[2], by=t.eff[i])
  th_sel <- th_new[i,indep_seq]
  print(ks.test(th_sel[1:as.integer(length(th_sel)/2)], th_sel[as.integer(length(th_sel)/2)+1:length(th_sel)]))
}



#Things to do:
#1. Change around the tuning parameters a bit to get 25% acceptance for the parameters 


#MCMC CORRELATION:
D <- rbind(th_new,B_new)

theta0 <- th[,n.mc];b0 <- b ## store starting values

V <- cov(t(D)) ## cov matrix of theta,b from first run
library(MASS) ## load mvrnorm MVN simulator



n.mc <- 50000 ## second sim length
pp <- mvrnorm(n.mc,rep(0,54),V) ## proposal perturbations

log_post_prop <- log_post <- log_posterior(theta0, y, b0, X, Z)
accept <- 0
th_corr <- matrix(0,4,n.mc); th_corr[,1] <- theta <- theta0
B_corr <- matrix(0,50,n.mc); b <- b0
for (i in 2:n.mc) { ## MH loop
  theta <- theta + pp[i,1:4]*0.3
  b <- b + pp[i,5:54]*0.01
  if (theta[4] > 0) log_post_prop <- log_posterior(theta, y, b, X, Z)
  if (theta[4] > 0 && runif(1) < exp(log_post_prop-log_post)) { 
    accept <- accept + 1
    log_post<-log_post_prop
  } else { ## reject
    theta <- th_corr[,i-1]
    b <- B_corr[,i-1]
  }
  th_corr[,i] <- theta;B_corr[,i] <- b
} ## end of MH loop
accept/n.mc


#Check for convergence and check acf...
par(mfrow=c(2,2))
for (j in 1:4) plot(th_corr[j,],type="l", xlab=variable_names[j])  
for (j in 1:4) acf(th_corr[j,])




#Produce plots to learn about the marginal densities of the parameters...
layout(matrix(c(1,2,1,2,3,4),2,3))

plot(1:n.mc,th_corr[1,],type="l",xlab="iteration",ylab="beta_0")
plot(1:n.mc,exp(th_corr[2,]),type="l",xlab="iteration",ylab="beta_1")
hist(th_corr[1,-(1:1000)],main="",xlab=expression(mu))
hist(exp(th_corr[2,-(1:1000)]),main="",xlab=expression(sigma)) 

plot(1:n.mc,th_corr[3,],type="l",xlab="iteration",ylab="beta_0")
plot(1:n.mc,exp(th_corr[4,]),type="l",xlab="iteration",ylab="beta_1")
hist(th_corr[3,-(1:1000)],main="",xlab=expression(mu))
hist(exp(th_corr[4,-(1:1000)]),main="",xlab=expression(sigma)) 


#Effective sample size, and checking for convergence using the Kolmogorov-Smirnov test. NOT WORKING - KS TEST SAYING NOT FROM SAME DISTRIBUTION.
n.eff <- c(0,0)
t.eff <- c(0,0)
par(mfrow=c(2,2))
for (j in 1:4) {
  acor <- acf(th_corr[j,])
  t.eff[j] <- 2*sum(acor[[1]]) + 1
  n.eff[j] <- dim(th_corr)[2]/t.eff[j]
}
n.eff<- as.integer(n.eff) #From Lec 9, this is the 'effective sample size' for the 4 parameters
t.eff<-as.integer(t.eff)


for (i in 1:4){
  indep_seq = seq(1, dim(th_corr)[2], by=t.eff[i])
  th_sel <- th_corr[i,indep_seq]
  print(ks.test(th_sel[1:as.integer(length(th_sel)/2)], th_sel[as.integer(length(th_sel)/2)+1:length(th_sel)]))
}



