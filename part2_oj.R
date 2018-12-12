
set.seed(528491)

fatigue <- read.table("https://people.bath.ac.uk/kai21/ASI/fatigue.txt")



N_prob_const_gamma <- deriv(
  expression(
    -log_sigma - log_alpha - delta*log(s - gamma) +
      ((1/exp(log_sigma)) - 1) * (log(y) - log_alpha - delta*log(s - gamma)) -
      (y / (exp(log_alpha) * (s - gamma)^delta))^(1/exp(log_sigma))
  ),
  namevec = c("log_alpha", "delta", "log_sigma"),
  function.arg = c("y", "s", "min_s", "gamma", "log_alpha", "delta", "log_sigma"),
  hessian = TRUE
)

N_surv_const_gamma <- deriv(
  expression(
    -(y / (exp(log_alpha) * (s - gamma)^delta))^(1/exp(log_sigma))
  ),
  namevec = c("log_alpha", "delta", "log_sigma"),
  function.arg = c("y", "s", "min_s", "gamma", "log_alpha", "delta", "log_sigma"),
  hessian = TRUE
)




# Negative log-likelihood
N_nll_const_gamma <- function(theta, y, gamma, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]

  -sum(
    (1-runout) * N_prob_const_gamma(y, stress, min(stress), gamma, log_alpha, delta, log_sigma),
    runout * N_surv_const_gamma(y, stress, min(stress), gamma, log_alpha, delta, log_sigma)
  )
}


# Gradient of negative log-likelihood
N_nll_gr_const_gamma <- function(theta, y, gamma, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]

  -colSums(rbind(
    (1-runout) * attr(N_prob_const_gamma(y, stress, min(stress), gamma, log_alpha, delta, log_sigma), "gradient"),
    runout * attr(N_surv_const_gamma(y, stress, min(stress), gamma, log_alpha, delta, log_sigma), "gradient")
  ))
}



# Initial params
theta0 <- c("log_alpha" = 2, "delta" = 1, "log_sigma" = 2)
gamma <- 70

# Minimise negative log-likelihood
N_opt_const_gamma <- optim(
  par = theta0,
  fn = N_nll_const_gamma,
  gr = N_nll_gr_const_gamma,
  y = fatigue$N,
  gamma = gamma,
  stress = fatigue$s,
  runout = fatigue$ro,
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 500)
)


# Calculate standard errors from inverse Hessian
N_std_err_const_gamma <- sqrt(diag(solve(N_opt$hessian)))

# 95% confidence interval for each parameter
# Note asymptotic distribution is normal
round(data.frame(
  val = N_opt_const_gamma$par,
  se = N_std_err_const_gamma,
  lower = N_opt_const_gamma$par - qnorm(0.975)*N_std_err_const_gamma,
  upper = N_opt_const_gamma$par + qnorm(0.975)*N_std_err_const_gamma
), 3)








# Expressions for Weibull probability density function and Weibull survival function
# l_gamma is logit (inverse sigmoid) transform of gamma, i.e.
#  l_gamma = log(gamma / (min(fatigue$s) - gamma))
N_prob <- deriv(
  expression(
      -log_sigma - log_alpha - delta*log(s - min_s/(1 + exp(-l_gamma))) +
        ((1/exp(log_sigma)) - 1) * (log(y) - log_alpha - delta*log(s - min_s/(1 + exp(-l_gamma)))) -
        (y / (exp(log_alpha) * (s - min_s/(1 + exp(-l_gamma)))^delta))^(1/exp(log_sigma))
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






log_posterior <- function(theta, y, b, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]
  mu_gamma <- theta[4]
  log_sigma_gamma <- theta[5]

  l_gamma <- log(b / (min(stress)-b))

  # Log conditional density of y given b
  lfy_b <- sum(
    (1-runout) * N_prob(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma),
    runout * N_surv(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma)
  )

  # Log marginal density of b
  lfb <- sum(dweibull(x = b, shape = 1/exp(log_sigma_gamma), scale = exp(mu_gamma), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb

  # Define log-prior (sum of log-priors for all parameters - here, just for log_sigma_b)
  log_prior <- dexp(exp(log_sigma_gamma), rate = 5, log = TRUE)

  # Log-posterior is log-prior plus log-likelihood
  lf + log_prior

}


log_posterior(c(1, 1, 1, 1, 1), fatigue$N, rep(0.1, 26), fatigue$s, fatigue$ro)




iters <- 100000
burnin <- 2000
pilot <- mcmc_mh(
  iters, burnin,
  c(log_alpha = 18, delta = -2, log_sigma = -1, mu_gamma = 4, log_sigma_gamma = -2),
  rep(66, 26), rep(0.03, 5), 0.12,
  fatigue$N, fatigue$s, fatigue$ro
)


D <- cbind(pilot$theta, pilot$b)[(burnin+1):iters, ]
cov_D <- cov(D)


adjusted <- mcmc_mh_cov(
  50000, 1000,
  drop(tail(pilot$theta, 1)),
  drop(tail(pilot$b, 1)), cov_D, 0.2,
  fatigue$N, fatigue$s, fatigue$ro
)



adjD <- cbind(adjusted$theta, adjusted$b)[-(1:burnin), ]
psych::pairs.panels(tail(adjD, 5000)[, 1:ncol(adjusted$theta)], pch = ".")


par(mfrow = c(2, 3))

ks_pvals <- vapply(colnames(adjusted$theta), function(nm) {

  # Calculate autocorrelation in MH sample for parameter
  y <- adjusted$theta[-(1:burnin), nm]
  autocorr <- acf(y, main = nm)

  # Autocorrelation length
  acl <- 2*sum(autocorr$acf) + 1

  # Sample of acl-spaced observations
  ac_smp <- y[seq(from = 1, to = length(y), by = acl)]

  # Test whether two halves of this sample are from same distribution
  # If p-value is significant, this means samples are (likely) from same dist
  idx <- sample(1:length(ac_smp), ceiling(length(ac_smp)/2), replace = FALSE)
  ks.test(ac_smp[idx], ac_smp[-idx])$p.value

}, FUN.VALUE = 0)

par(mfrow = c(1, 1))

# p-values as just calculated
ks_pvals



# Credible intervals for each parameter
cred_ints <- vapply(colnames(adjusted$theta), function(nm) {

  y <- adjusted$theta[-(1:burnin), nm]
  quantile(y, c(0.025, 0.975))

}, FUN.VALUE = c(0, 0))

cred_ints











intgrl <- function(theta) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]
  mu_gamma <- theta[4]
  log_sigma_gamma <- theta[5]

  vapply(seq_len(nrow(fatigue)), function(k) {
    integrate(
      function(gamma) {
        ((1-fatigue$ro[k]) * dweibull(fatigue$N[k], shape = 1/exp(log_sigma), scale = (exp(log_alpha)*(fatigue$s[k] - gamma)^delta)) +
           fatigue$ro[k] * pweibull(fatigue$N[k], shape = 1/exp(log_sigma), scale = (exp(log_alpha)*(fatigue$s[k] - gamma)^delta), lower.tail = FALSE)) *
          dweibull(gamma, shape = 1/exp(log_sigma_gamma), scale = exp(mu_gamma))
      },
      lower = 0, upper = fatigue$s[k]
    )$value
  }, FUN.VALUE = 0)

}


intgrl(c(18, -2, -1, 4, -2))


# Calculate log posterior again but now without any priors

log_post_no_prior <- function(theta, y, b, stress, runout) {

  log_alpha <- theta[1]
  delta <- theta[2]
  log_sigma <- theta[3]
  mu_gamma <- theta[4]
  log_sigma_gamma <- theta[5]

  l_gamma <- log(b / (min(stress)-b))

  # Log conditional density of y given b
  lfy_b <- sum(
    (1-runout) * N_prob(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma),
    runout * N_surv(y, stress, min(stress), log_alpha, delta, log_sigma, l_gamma)
  )

  # Log marginal density of b
  lfb <- sum(dweibull(x = b, shape = 1/exp(log_sigma_gamma), scale = exp(mu_gamma), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lfy_b + lfb
}


# Create full proposal distribution with denominator included
log_posterior <- function(theta, y, b, stress, runout) {
  l_fn <- sum(log(intgrl(theta)))
  l_fn_b_fb <- log_post_no_prior(theta, y, b, stress, runout)
  l_fn_b_fb - l_fn
}

#MCMC again (still taking correlation into account) with full posterior
adjusted_new <- mcmc_mh_cov(
  50000, 1000,
  drop(tail(pilot$theta, 1)),
  drop(tail(pilot$b, 1)), cov_D, 0.25,
  fatigue$N, fatigue$s, fatigue$ro
)
