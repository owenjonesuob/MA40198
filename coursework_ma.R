rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

nll <- function(theta){

  l1 = dweibull(rats$time[rats$status==1&rats$rx==1], shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]), log=TRUE)
  l2 = dweibull(rats$time[rats$status==1&rats$rx==0], shape=1/exp(theta[3]), scale=exp(theta[1]), log=TRUE)

  l3 = pweibull(rats$time[rats$status==0&rats$rx==1], shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]), log.p=TRUE, lower.tail=FALSE)
  l4 = pweibull(rats$time[rats$status==0&rats$rx==0], shape=1/exp(theta[3]), scale=exp(theta[1]), log.p=TRUE, lower.tail=FALSE)

  -(sum(l1)+sum(l2)+sum(l3)+sum(l4))

}


theta0=c(1,1,1)

q=optim(par=theta0,fn=nll, method="BFGS", hessian=TRUE)
inv_hess = solve(q$hessian)

#2.
q$par[2] + 1.96*sqrt(inv_hess[2,2])
q$par[2] - 1.96*sqrt(inv_hess[2,2])

#Interpretation: This shows that at the 5% significance level, beta1 is significantly different from 0 and therefore
#we can say that there is sufficient evidence to suggest that recieving treatment has an effect on the time it takes
#for a tumour to appear. We have a negative estimate, which decreases the scale parameter in the Weibull distribution.
#This in turn DECREASES/INCREASES the esimated times for tumour appearance.

#3.

install.packages("actuar")
library(actuar)

nll_ll <- function(theta){

  shape = 1/exp(theta[3])
  scale_0 = exp(theta[1])
  scale_1 = exp(theta[1]+theta[2])

  l1 = log(shape)-log(scale_1)+(shape-1)*log(rats$time[rats$status==1&rats$rx==1]/scale_1) - 2*log(1+(rats$time[rats$status==1&rats$rx==1]/scale_1)^shape)
  l2 = log(shape)-log(scale_0)+(shape-1)*log(rats$time[rats$status==1&rats$rx==0]/scale_0) - 2*log(1+(rats$time[rats$status==1&rats$rx==0]/scale_0)^shape)

  l3 = -log(1+(rats$time[rats$status==0&rats$rx==1]/scale_1)^shape)
  l4 = -log(1+(rats$time[rats$status==0&rats$rx==0]/scale_0)^shape)

  -(sum(l1)+sum(l2)+sum(l3)+sum(l4))

}


theta0_ll=c(1,1,1)
q_ll=optim(par=theta0_ll,fn=nll_ll, method="BFGS", hessian=TRUE)
inv_hess_ll = solve(q_ll$hessian)
q_ll$par[2] + 1.96*sqrt(inv_hess_ll[2,2])
q_ll$par[2] - 1.96*sqrt(inv_hess_ll[2,2])


#4.
rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1,rats)

X <- model.matrix(~ rx, rats)


lfyb <- function(theta, y, b, X, Z){
  eta_0 = as.numeric(theta[1] + Z%*%b)
  eta_1 = as.numeric(X%*%theta[1:2] + Z%*%b)
  scale_0 = exp(eta_0)
  scale_1 = exp(eta_1)
  shape = 1/exp(theta[3])

  shape = 1/exp(theta[3])

  l1 = dweibull(y[rats$status==1&rats$rx==1], shape=shape, scale=scale_1, log=TRUE)
  l2 = dweibull(y[rats$status==1&rats$rx==0], shape=shape, scale=scale_0, log=TRUE)

  l3 = pweibull(y[rats$status==0&rats$rx==1], shape=shape, scale=scale_1, log.p=TRUE, lower.tail=FALSE)
  l4 = pweibull(y[rats$status==0&rats$rx==0], shape=shape, scale=scale_0, log.p=TRUE, lower.tail=FALSE)

  lfy_b = sum(l1)+sum(l2)+sum(l3)+sum(l4)
  lfb = sum(dnorm(b, 0, exp(theta[4]), log=TRUE))
  -lfy_b - lfb


}

lal <- function(theta, y, X, Z){

  b <- rep(0,ncol(Z))
  opt = optim(par=b, lfyb, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)


  la <- -lfyb(theta, y, opt$par, X, Z) + length(b)*log(2*pi)/2 - sum(log(abs(diag(solve(opt$hessian)))))/2

  attr(la,"b") <- as.numeric(b)

  -la

}

theta0=c(1.1749184,  0.6614587, -1.4266831,  1.7078852)

optim(par=theta0, lal, X=X, y=rats$time, Z=Z, method="BFGS", hessian=TRUE)

#works w/ nelder mead



#5

#Posterior proportional to prior x likelihood.
#Doing log posterior for ease of calculations
posterior <- function(theta, y, b, X, Z){
  eta_0_0 = as.numeric(theta[1] + Z%*%b)[rats$status==0&rats$rx==0]
  # treatment==1 means include beta1 in eta
  # rx==1 means use survival function rather than pdf
  eta_1_0 = as.numeric(X%*%theta[1:2] + Z%*%b)[rats$status==0&rats$rx==1]

  eta_0_1 = as.numeric(theta[1] + Z%*%b)[rats$status==1&rats$rx==0]
  eta_1_1 = as.numeric(X%*%theta[1:2] + Z%*%b)[rats$status==1&rats$rx==1]

  scale_0_0 = exp(eta_0_0)
  scale_0_1 = exp(eta_0_1)

  scale_1_0 = exp(eta_1_0)
  scale_1_1 = exp(eta_1_1)

  shape = 1/exp(theta[3])


  l1 = dweibull(y[rats$status==1&rats$rx==1], shape=shape, scale=scale_1_1, log=TRUE)
  l2 = dweibull(y[rats$status==0&rats$rx==1], shape=shape, scale=scale_1_0, log=TRUE)

  l3 = pweibull(y[rats$status==1&rats$rx==0], shape=shape, scale=scale_0_1, log.p=TRUE, lower.tail=FALSE)
  l4 = pweibull(y[rats$status==0&rats$rx==0], shape=shape, scale=scale_0_0, log.p=TRUE, lower.tail=FALSE)

  lfy_b = sum(l1)+sum(l2)+sum(l3)+sum(l4)
  lfb = sum(dnorm(b, 0, exp(theta[4]), log=TRUE))
  likelihood = lfy_b + lfb
  sigma_prior = dexp(exp(theta[4]), rate=5, log=TRUE)

  likelihood + sigma_prior
  #No other priors are needed here as they are all assumed to be improper uniform, i.e. proportional to 1

}


#OWEN'S WORK:
#

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




posterior <- function(theta, y, b, X, Z) {

  beta <- theta[1:2]
  log_sigma <- theta[3]
  log_sigma_b <- theta[4]

  eta <- X%*%beta + Z%*%b

  status <- X[, 2]

  # Log conditional density of y given b
  # Sometimes we sum lots of -Infs and get NaN... prevent this
  lfy_b <- max(sum(
    status * weib_re_prob(y, eta, log_sigma),
    (1-status) * weib_re_surv(y, eta, log_sigma)
  ), -Inf, na.rm = TRUE)

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb

  # Define prior (sum of priors for all parameters - here, just for log_sigma_b)
  prior <- dexp(exp(log_sigma_b), rate = 5, log = TRUE)

  # (Log) posterior is (log) prior plus (log) likelihood
  lf + prior

}


posterior(c(1, 1, 1, 1), rats$time, rep(0.1, ncol(Z)), X, Z)


#QQQ ISSUE AT THE MOMENT - NAN'S ARE PRODUCED........




mcmc_mh <- function(iters, burnin, init_params, tuners, y, X, Z, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = 4)
  theta_vals[1, ] <- init_params

  acceptance <- rep(NA, iters)

  pb <- txtProgressBar(min = 0, max = iters, style = 3)

  for (k in seq_len(iters)) {

    theta_curr <- theta_vals[k, ]
    theta_prop <- rnorm(4, mean = theta_curr, sd = tuners)

    # Random effects have standard deviation exp(log_sigma_b)
    b_curr <- rnorm(ncol(Z), mean = 0, sd = exp(theta_curr[4]))
    b_prop <- rnorm(ncol(Z), mean = 0, sd = exp(theta_prop[4]))


    post_curr <- posterior(theta_curr, y, b_curr, X, Z)
    post_prop <- posterior(theta_prop, y, b_prop, X, Z)


    accept_prob <- exp(post_prop - post_curr)

    if (accept_prob > runif(1)) {

      theta_vals[k+1, ] <- theta_prop
      acceptance[k] <- 1

    } else {

      theta_vals[k+1, ] <- theta_curr
      acceptance[k] <- 0
    }

    setTxtProgressBar(pb, value = k, label = if (k <= burnin) "Burn-in" else "Main run")

  }


  if (show_plot) {
    par(mfrow = c(2, 2))
    apply(theta_vals, 2, function(x) {
      plot(x, type = "l")
      abline(v = burnin, col = "red")
    })
    par(mfrow = c(1, 1))
  }

  close(pb)
  cat(sprintf("Total acceptance:     %.3f%%\n", mean(acceptance)*100))
  cat(sprintf("Burned-in acceptance: %.3f%%\n", mean(acceptance[-(1:burnin)])*100))

  theta_vals
}



zz <- mcmc_mh(5000, 3000, rep(1, 4), c(0.5, 0.5, 0.5, 0.1), rats$time, X, Z)









for (tuning_param in seq(from = 1, to = 2, by = 0.1)) {

  # Initial guess
  theta_init <- c(0.5, 0.5, 2, 2)


  mcmc_mh(10000, 0, tuning_param, theta_init, rats$time, X, Z, show_plot = FALSE)


  print(tuning_param)
  print(mean(accepted))
}



#Use the best tuning parameter from above and run a full MCMC of length 100000

n=100000
tuning_param= 0.016

theta_vals = matrix(0,n+1,4)
theta_vals[1,]=c(0.5,0.5,2, 2)
accepted = c()
for (i in 1:n){
  theta_proposed = c(0,0,0,0)
  theta_proposed[1] =  rnorm(1,mean=theta_vals[i,1],sd=tuning_param) #beta_0
  theta_proposed[2] =  rnorm(1,mean=theta_vals[i,2],sd=tuning_param) #beta_1
  theta_proposed[3] =  rnorm(1,mean=theta_vals[i,3],sd=tuning_param) #log(sig)
  theta_proposed[4] =  rnorm(1,mean=theta_vals[i,4],sd=tuning_param) #log(sig.b)

  theta_new = c(theta_proposed[1],theta_proposed[2],exp(theta_proposed[3]), exp(theta_proposed[4]))
  theta_old = c(theta_vals[i,1],theta_vals[i,1],exp(theta_vals[i,3]), exp(theta_vals[i,4]))

  b_new = rnorm(50, mean=0, sd=theta_new[4])
  b_old = rnorm(50, mean=0, sd=theta_old[4])

  log_new_val = posterior(theta_new, y, b_new, X, Z)
  log_old_val = posterior(theta_old, y, b_old, X, Z)

  accept_prob = min(0,log_new_val-log_old_val)
  alpha <- log(runif(1))
  if (accept_prob >= alpha){
    accepted[i]=1
    theta_vals[i+1,] = c(theta_proposed[1],theta_proposed[2],theta_proposed[3], theta_proposed[4])
  } else {
    accepted[i] = 0
    theta_vals[i+1,] = theta_vals[i,]
  }
}
print(sum(accepted)/length(accepted))

#Plots of values to check for convergence
par(mfrow=c(2,2))
plot(theta_vals[,1])
plot(theta_vals[,2])
plot(theta_vals[,3])
plot(theta_vals[,4])





#NEXT PART: USING THE COVARIANCE BETWEEN THE PARAMETERS... NEED TO WORK ON THIS!
library(mvtnorm)


burned_in = theta_vals[50000:dim(theta_vals)[1],]
#ESTIMATED COVARIANCE MATRIX:
sig.est = cov(burned_in)
mu.est = apply(burned_in, 2, mean)
#Use p158 v2 'shrunken version'

n=10000
tuning_param= 0.016

theta_vals = matrix(0,n+1,4)
theta_vals[1,]=c(0.5,0.5,2, 2)
accepted = c()
for (i in 1:n){
  theta_proposed = rmvnorm(4, mean=mu.est, sigma=sig.est)

  theta_new = c(theta_proposed[1],theta_proposed[2],exp(theta_proposed[3]), exp(theta_proposed[4]))
  theta_old = c(theta_vals[i,1],theta_vals[i,1],exp(theta_vals[i,3]), exp(theta_vals[i,4]))

  b_new = rnorm(50, mean=0, sd=theta_new[4])
  b_old = rnorm(50, mean=0, sd=theta_old[4])

  log_new_val = posterior(theta_new, y, b_new, X, Z)
  log_old_val = posterior(theta_old, y, b_old, X, Z)

  accept_prob = min(0,log_new_val-log_old_val+dmvnorm(theta_old, mean=mu.est, sigma=sig.est, log=TRUE) - dmvnorm(theta_new, mean=mu.est, sigma=sig.est, log=TRUE)) #***NOW NEED TO INCLUDE PROPOSAL DISTRIBUTIONS***
  alpha <- log(runif(1))
  if (accept_prob >= alpha){
    accepted[i]=1
    theta_vals[i+1,] = c(theta_proposed[1],theta_proposed[2],theta_proposed[3], theta_proposed[4])
  } else {
    accepted[i] = 0
    theta_vals[i+1,] = theta_vals[i,]
  }
}

par(mfrow=c(2,2))
plot(theta_vals[,1])
plot(theta_vals[,2])
plot(theta_vals[,3])
plot(theta_vals[,4])


#STILL TO DO:
#FIGURE OUT WHAT TO DO ABOUT B AND JUSTIFY IT
#PRODUCE PLOTS AND JUSTIFICATIONS FOR CONVERGENCE
#THINK OF GOOD STARTING POINTS AND JUSTIFY THEM
#JUSTIFY SOME GOOD TUNING PARAMETERS FOR THE 4 PARAMETERS
#Use marginal likelihoods as starting values and use their SDs as the tuning parameters....
#FIGURE OUT HOW TO DO THE COVARIANCE BETWEEN PARAMETERS STUFF - GETTING THERE...
#'PRODUCE PLOTS TO LEARN ABOUT THE SHAPE OF THE MARGINAL POSTERIOR DENSITIES OF THE PARAMETERS'??









#SECTION 2
#QUESTION 1

fatigue <- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")

gamma_try = min(fatigue$s[fatigue$ro==0])-1

nll_s2 <- function(theta){
  sigma = exp(theta[3])
  alpha = exp(theta[1])
  delta = theta[2]
  not_censored = sum(dweibull(fatigue$N[fatigue$ro==0], scale = (alpha*(fatigue$s[fatigue$ro==0] - gamma_try)^delta), shape=(1/sigma), log=TRUE))
  censored = sum(pweibull(fatigue$N[fatigue$ro==1], scale = (alpha*(fatigue$s[fatigue$ro==1] - gamma_try)^delta), shape=(1/sigma), log.p=TRUE, lower.tail=FALSE))
  -not_censored-censored
}

theta0 = c(1,1,1)
optim(theta0, nll_s2, method="BFGS", hessian=TRUE, control=list(trace=1, REPORT=1))

q=optim(theta0, nll_s2, control=list(trace=1, REPORT=1))
q$par


#Gamma between 0 and 80
#2

nll_s2_2 <- function(theta){
  gamma = 80*(1/(1+exp(-theta[4])))
  sigma = exp(theta[3])
  alpha = exp(theta[1])
  delta = theta[2]
  not_censored = sum(dweibull(fatigue$N[fatigue$ro==0], scale = (alpha*(fatigue$s[fatigue$ro==0] - gamma_try)^delta), shape=(1/sigma), log=TRUE))
  censored = sum(pweibull(fatigue$N[fatigue$ro==1], scale = (alpha*(fatigue$s[fatigue$ro==1] - gamma_try)^delta), shape=(1/sigma), log.p=TRUE, lower.tail=FALSE))
  -not_censored-censored

}

theta0 = c(10,10,10,10)
optim(theta0, nll_s2_2,control=list(trace=1, REPORT=1),method="BFGS", hessian=TRUE)
