
# Import data
rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")



# Negative log-likelihood using values from built-in weibull functions
weib_nll_numeric <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]


  shape <- 1/exp(log_sigma)
  scale <- exp(beta0 + beta1*treated)

  -sum(
    status * dweibull(t, shape = shape, scale = scale, log = TRUE),
    (1-status) * pweibull(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  )
}



# Expressions for Weibull probability density function and Weibull survival function
weib_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/exp(beta0 + beta1*treated) * (t/exp(beta0 + beta1*treated))^((1/exp(log_sigma))-1) * exp(-(t/exp(beta0 + beta1*treated))^(1/exp(log_sigma))))
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)


weib_surv <- deriv(
  expression(
    -(t/exp(beta0 + beta1*treated))^(1/exp(log_sigma))
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)



# Negative log-likelihood
weib_nll <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -sum(
    status * weib_prob(t, beta0, beta1, log_sigma, treated),
    (1-status) * weib_surv(t, beta0, beta1, log_sigma, treated)
  )
}


# Gradient of negative log-likelihood
weib_nll_gr <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -colSums(rbind(
    status * attr(weib_prob(t, beta0, beta1, log_sigma, treated), "gradient"),
    (1-status) * attr(weib_surv(t, beta0, beta1, log_sigma, treated), "gradient")
  ))
}



# Initial params
theta0 <- c("beta0" = 2, "beta1" = 2, "log_sigma" = 2)


# Minimise negative log-likelihood
weib_opt <- optim(
  par = theta0,
  fn = weib_nll,
  gr = weib_nll_gr,
  t = rats$time,
  treated = rats$rx,
  status = rats$status,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)


# Calculate standard errors from inverse Hessian
weib_std_err <- sqrt(diag(solve(weib_opt$hessian)))
names(weib_std_err) <- c("beta0", "beta1", "log_sigma")


# 95% confidence interval for beta1
# Note asymptotic distribution is normal (for MLE estimate)
weib_opt$par[2] + c(-1, 1)*qnorm(0.975)*weib_std_err[2]









# Expressions for log-logistic probability density function and log-logistic survival function
llog_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/exp(beta0 + beta1*treated) * (t/exp(beta0 + beta1*treated))^((1/exp(log_sigma))-1) / (1 + (t/exp(beta0 + beta1*treated))^(1/exp(log_sigma)))^2)
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)


llog_surv <- deriv(
  expression(
    -log(1 + (t/exp(beta0 + beta1*treated))^(1/exp(log_sigma)))
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)



# Negative log-likelihood
llog_nll <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -sum(
    status * llog_prob(t, beta0, beta1, log_sigma, treated),
    (1-status) * llog_surv(t, beta0, beta1, log_sigma, treated)
  )
}


# Gradient of negative log-likelihood
llog_nll_gr <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -colSums(rbind(
    status * attr(llog_prob(t, beta0, beta1, log_sigma, treated), "gradient"),
    (1-status) * attr(llog_surv(t, beta0, beta1, log_sigma, treated), "gradient")
  ))
}




# Minimise negative log-likelihood
llog_opt <- optim(
  par = theta0,
  fn = llog_nll,
  gr = llog_nll_gr,
  t = rats$time,
  treated = rats$rx,
  status = rats$status,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)


# Calculate standard errors from inverse Hessian
llog_std_err <- sqrt(diag(solve(llog_opt$hessian)))

# 95% confidence interval for each parameter
# Note asymptotic distribution is normal
data.frame(par = names(llog_opt$par),
           val = llog_opt$par,
           se = llog_std_err,
           lower = llog_opt$par - qnorm(0.975)*llog_std_err,
           upper = llog_opt$par + qnorm(0.975)*llog_std_err
)







X <- model.matrix(~ 1 + rx, data = rats)
Z <- model.matrix(~ litter - 1, data = rats)




weib_me_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/exp(eta) * (t/exp(eta))^((1/exp(log_sigma))-1) * exp(-(t/exp(eta))^(1/exp(log_sigma))))
  ),
  namevec = c("eta", "log_sigma"),
  function.arg = c("t", "b", "eta", "log_sigma"),
  hessian = TRUE
)


weib_me_surv <- deriv(
  expression(
    -(t/exp(eta))^(1/exp(log_sigma))
  ),
  namevec = c("eta", "log_sigma"),
  function.arg = c("t", "b", "eta", "log_sigma"),
  hessian = TRUE
)




lfyb <- function(theta, y, b, X, Z) {

  beta <- theta[1:2]
  log_sigma <- theta[3]
  log_sigma_b <- theta[4]

  eta <- X%*%beta + Z%*%b

  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_me_prob(t, b, eta, log_sigma),
    (1-status) * weib_me_surv(t, b, eta, log_sigma)
  )

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb

  # Now gradient and Hessian
  g <- colSums(rbind(
    status * attr(weib_me_prob(t, eta, log_sigma, treated), "gradient"),
    (1-status) * attr(weib_me_surv(t, eta, log_sigma, treated), "gradient")
  ))

  H <- colSums(rbind(
    status * attr(weib_me_prob(t, eta, log_sigma, treated), "hessian"),
    (1-status) * attr(weib_me_surv(t, eta, log_sigma, treated), "hessian")
  ))



}



