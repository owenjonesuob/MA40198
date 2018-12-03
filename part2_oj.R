
fatigue <- read.table("https://people.bath.ac.uk/kai21/ASI/fatigue.txt")


# Expressions for Weibull probability density function and Weibull survival function
# l_gamma is logit (inverse sigmoid) transform of gamma, i.e.
#  l_gamma = log(gamma / (min(fatigue$s) - gamma))
N_prob <- deriv(
  expression(
    log(
      (1/exp(log_sigma)) / (exp(log_alpha) * (s - min_s/(1 + exp(-l_gamma)))^delta) *
        (y / (exp(log_alpha) * (s - min_s/(1 + exp(-l_gamma)))^delta))^((1/exp(log_sigma)) - 1) *
        exp(-(y / (exp(log_alpha) * (s - min_s/(1 + exp(-l_gamma)))^delta))^(1/exp(log_sigma)))
    )
  ),
  namevec = c("log_alpha", "delta", "log_sigma", "l_gamma"),
  function.arg = c("y", "s", "min_s", "log_alpha", "delta", "log_sigma", "l_gamma"),
  hessian = TRUE
)


N_surv <- deriv(
  expression(
    -(y / (exp(log_alpha) * (s - min_s/(1 + exp(-l_gamma)))^delta))^(1/exp(log_sigma))
  ),
  namevec = c("log_alpha", "delta", "log_sigma", "l_gamma"),
  function.arg = c("y", "s", "min_s", "log_alpha", "delta", "log_sigma", "l_gamma"),
  hessian = TRUE
)



# Negative log-likelihood
N_nll <- function(theta, y, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]
  l_gamma <- theta[4]

  -sum(
    (1-runout) * N_prob(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma),
    runout * N_surv(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma)
  )
}


# Gradient of negative log-likelihood
N_nll_gr <- function(theta, y, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]
  l_gamma <- theta[4]

  -colSums(rbind(
    (1-runout) * attr(N_prob(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma), "gradient"),
    runout * attr(N_surv(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma), "gradient")
  ))
}



# Initial params
# Certain combinations seem to lead to strange (non-optimal) results...
# e.g. c(1, 1, 1, 1)
theta0 <- c("log_alpha" = 2, "delta" = 1, "log_sigma" = 2, "l_gamma" = 0)

# Minimise negative log-likelihood
N_opt <- optim(
  par = theta0,
  fn = N_nll,
  gr = N_nll_gr,
  y = fatigue$N,
  stress = fatigue$s,
  runout = fatigue$ro,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)


# Calculate standard errors from inverse Hessian
N_std_err <- sqrt(diag(solve(N_opt$hessian)))

# 95% confidence interval for each parameter
# Note asymptotic distribution is normal
round(data.frame(
  val = N_opt$par,
  se = N_std_err,
  lower = N_opt$par - qnorm(0.975)*N_std_err,
  upper = N_opt$par + qnorm(0.975)*N_std_err
), 3)


# Recover gamma
min(fatigue$s) / (1 + exp(-N_opt$par[4]))






alpha <- exp(N_opt$par[1])
delta <- N_opt$par[2]
sigma <- exp(N_opt$par[3])
gamma <- min(fatigue$s) / (1 + exp(-N_opt$par[4]))


library(ggplot2)

ggplot(fatigue, aes(x = s)) +
  geom_point(aes(y = N)) +
  stat_function(fun = function(s) {alpha * (s - gamma)^delta * qweibull(0.1, shape = 1, scale = 1)}, aes(color = "10% quantile")) +
  stat_function(fun = function(s) {alpha * (s - gamma)^delta * qweibull(0.5, shape = 1, scale = 1)}, aes(color = "50% quantile"), linetype = "dashed") +
  ggtitle("Estimated quantiles of fatigue") + labs(x = "Stress", y = "N", color = "Estimated quantiles") +
  theme_bw() + theme(legend.position = c(0.8, 0.85), legend.background = element_rect(color = "black"))
