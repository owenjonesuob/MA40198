
fatigue <- read.table("https://people.bath.ac.uk/kai21/ASI/fatigue.txt")


# Expressions for Weibull probability density function and Weibull survival function
N_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/(exp(log_alpha)*(s - gamma)^delta) * (y/(exp(log_alpha)*(s - gamma)^delta))^((1/exp(log_sigma))-1) * exp(-(y/(exp(log_alpha)*(s - gamma)^delta))^(1/exp(log_sigma))))
  ),
  namevec = c("log_alpha", "delta", "log_sigma"),
  function.arg = c("y", "s", "gamma", "log_alpha", "delta", "log_sigma"),
  hessian = TRUE
)


N_surv <- deriv(
  expression(
    -(y/(exp(log_alpha)*(s - gamma)^delta))^(1/exp(log_sigma))
  ),
  namevec = c("log_alpha", "delta", "log_sigma"),
  function.arg = c("y", "s", "gamma", "log_alpha", "delta", "log_sigma"),
  hessian = TRUE
)



# Negative log-likelihood
N_nll <- function(theta, y, stress, runout, gamma) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]

  -sum(
    (1-runout) * N_prob(y, stress, gamma, log_alpha, delta, log_sigma),
    runout * N_surv(y, stress, gamma, log_alpha, delta, log_sigma)
  )
}


# Gradient of negative log-likelihood
N_nll_gr <- function(theta, y, stress, runout, gamma) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]

  -colSums(rbind(
    (1-runout) * attr(N_prob(y, stress, gamma, log_alpha, delta, log_sigma), "gradient"),
    runout * attr(N_surv(y, stress, gamma, log_alpha, delta, log_sigma), "gradient")
  ))
}



# Initial params
theta0 <- c("log_alpha" = 2, "delta" = -2, "log_sigma" = 2)
gamma <- 79.3#runif(1, min = 0, max = min(fatigue$s))

# Minimise negative log-likelihood
N_opt <- optim(
  par = theta0,
  fn = N_nll,
  gr = N_nll_gr,
  y = fatigue$N,
  stress = fatigue$s,
  runout = fatigue$ro,
  gamma = gamma,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)


# Calculate standard errors from inverse Hessian
N_std_err <- sqrt(diag(solve(N_opt$hessian)))

# 95% confidence interval for each parameter
# Note asymptotic distribution is normal
data.frame(
  par = names(N_opt$par),
  val = N_opt$par,
  se = N_std_err,
  lower = N_opt$par - qnorm(0.975)*N_std_err,
  upper = N_opt$par + qnorm(0.975)*N_std_err
)
