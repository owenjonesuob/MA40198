
# Import data
rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")



# Negative log-likelihood using values from built-in weibull functions
nll_numeric <- function(theta, t, treated, status) {

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
d_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/exp(beta0 + beta1*treated) * (t/exp(beta0 + beta1*treated))^((1/exp(log_sigma))-1) * exp(-(t/exp(beta0 + beta1*treated))^(1/exp(log_sigma))))
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)


d_surv <- deriv(
  expression(
    -(t/exp(beta0 + beta1*treated))^(1/exp(log_sigma))
  ),
  namevec = c("beta0", "beta1", "log_sigma"),
  function.arg = c("t", "beta0", "beta1", "log_sigma", "treated"),
  hessian = TRUE
)



# Negative log-likelihood
nll <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -sum(
    status * d_prob(t, beta0, beta1, log_sigma, treated),
    (1-status) * d_surv(t, beta0, beta1, log_sigma, treated)
  )
}


# Gradient of negative log-likelihood
nll_gr <- function(theta, t, treated, status) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]

  -colSums(rbind(
    status * attr(d_prob(t, beta0, beta1, log_sigma, treated), "gradient"),
    (1-status) * attr(d_surv(t, beta0, beta1, log_sigma, treated), "gradient")
  ))
}



# Initial params
theta0 <- c("beta0" = 2, "beta1" = 2, "log_sigma" = 2)


# Minimise negative log-likelihood
opt <- optim(
  par = theta0,
  fn = nll,
  gr = nll_gr,
  t = rats$time,
  treated = rats$rx,
  status = rats$status,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)


# Calculate standard errors from inverse Hessian
std_err <- sqrt(diag(solve(opt$hessian)))
names(std_err) <- c("beta0", "beta1", "log_sigma")


# 95% confidence interval for beta1
# Note asymptotic distribution is
opt$par[2] + c(-1, 1)*qnorm(0.975)*std_err[2]
