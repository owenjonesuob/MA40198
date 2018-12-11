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
  log_sigma_b <- theta[4]

  eta <- X%*%beta + Z%*%b

  time <- y[, "time"]
  status <- y[, "status"]

  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_re_prob(time, eta, log_sigma),
    (1-status) * weib_re_surv(time, eta, log_sigma)
  )

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb

  # Define log-prior (sum of log-priors for all parameters - here, just for log_sigma_b)
  log_prior <- dexp(exp(log_sigma_b), rate = 5, log = TRUE)

  # Log-posterior is log-prior plus log-likelihood
  lf + log_prior

}


log_posterior(c(1, 1, 1, 1), rats[, c("time", "status")], rep(0.1, ncol(Z)), X, Z)





mcmc_mh <- function(iters, burnin, init_params, tuners, b_tuner, y, X, Z, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = 4)
  theta_vals[1, ] <- init_params

  b_vals <- matrix(NA, nrow = iters+1, ncol = ncol(Z))
  b_vals[1, ] <- 0

  acceptance <- rep(NA, iters)
  acceptance_b <- rep(NA, iters)

  log_post <- log_posterior(init_params, y, b_vals[1, ], X, Z)

  # Progress bar in console
  pb <- txtProgressBar(min = 0, max = iters, style = 3)


  # MCMC loop
  for (k in seq_len(iters)) {


    # "Non-random" parameters
    theta_prop <- rnorm(4, mean = theta_vals[k, ], sd = tuners)


    log_post_prop <- log_posterior(theta_prop, y, b_vals[k, ], X, Z)

    accept_prob <- exp(log_post_prop - log_post)

    if (accept_prob > runif(1)) {

      theta_vals[k+1, ] <- theta_prop
      log_post <- log_post_prop
      acceptance[k] <- TRUE

    } else {

      theta_vals[k+1, ] <- theta_vals[k, ]
      acceptance[k] <- FALSE
    }


    # "Random" parameters (proposed/accepted separately)

    # Now hold other parameters steady
    # Random effects have standard deviation exp(log_sigma_b)
    b_prop <- rnorm(ncol(Z), mean = b_vals[k, ], sd = b_tuner)

    # Use (possibly newly accepted) values of theta
    log_post_prop_b <- log_posterior(theta_vals[k+1, ], y, b_prop, X, Z)

    accept_prob_b <- exp(log_post_prop_b - log_post)

    if (accept_prob_b > runif(1)) {

      b_vals[k+1, ] <- b_prop
      log_post <- log_post_prop_b
      acceptance_b[k] <- TRUE

    } else {

      b_vals[k+1, ] <- b_vals[k, ]
      acceptance_b[k] <- FALSE
    }

    setTxtProgressBar(pb, value = k)
  }


  if (show_plot) {
    par(mfrow = c(2, 2))
    apply(theta_vals, 2, function(x) {
      plot(x, type = "l")
      abline(v = burnin, col = "red")
      rect(0, min(x), burnin, min(y), density = 10, col = "red")
    })
    par(mfrow = c(1, 1))
  }

  close(pb)
  cat(sprintf("Total acceptance:       %2.3f%%\n", mean(acceptance)*100))
  cat(sprintf("Burned-in acceptance:   %2.3f%%\n", mean(acceptance[-(1:burnin)])*100))
  cat(sprintf("Burned-in b acceptance: %2.3f%%\n", mean(acceptance_b[-(1:burnin)])*100))

  list(theta = theta_vals, b = b_vals)
}


zz <- mcmc_mh(50000, 2000, c(4, 0, 0, -1), c(0.1, 0.1, 0.1, 0.1), 0.07, rats[, c("time", "status")], X, Z)



pairs(cbind(zz$theta, zz$b[, 1:4]), pch = ".")


# TODO below this




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
