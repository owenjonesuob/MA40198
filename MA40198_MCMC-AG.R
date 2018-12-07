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
for (j in 1:3) plot(th[j,],type="l") 

