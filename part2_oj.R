
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


log_posterior(c(1, 1, 1, 1, 1), fatigue$N, rep(0.1, 26), fatigue$s, fatigue$runout)







mcmc_mh_f <- function(iters, burnin, init_params, init_bs, tuners, b_tuner, y, stress, runout, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = 5)
  theta_vals[1, ] <- init_params

  b_vals <- matrix(NA, nrow = iters+1, ncol = 26)
  b_vals[1, ] <- init_bs

  acceptance <- rep(NA, iters)
  acceptance_b <- rep(NA, iters)

  log_post <- log_posterior(init_params, y, b_vals[1, ], stress, runout)

  # Progress bar in console
  pb <- txtProgressBar(min = 0, max = iters, style = 3)


  # MCMC loop
  for (k in seq_len(iters)) {


    # "Non-random" parameters
    theta_prop <- rnorm(5, mean = theta_vals[k, ], sd = tuners)


    log_post_prop <- log_posterior(theta_prop, y, b_vals[k, ], stress, runout)

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
    b_prop <- rnorm(26, mean = b_vals[k, ], sd = b_tuner)

    # Use (possibly newly accepted) values of theta
    log_post_prop_b <- log_posterior(theta_vals[k+1, ], y, b_prop, stress, runout)

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
    par(mfrow = c(2, 3))
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


pilot <- mcmc_mh_f(50000, 15000, c(18, -2, -1, 4, -2), rep(66, 26), rep(0.025, 5), 0.2, fatigue$N, fatigue$s, fatigue$ro, FALSE)


D <- cbind(pilot$theta, pilot$b)[15001:50001, ]

cov_D <- cov(D)






mcmc_mh_cov_f <- function(iters, burnin, init_params, init_bs, cov_matrix, tuner, y, stress, runout, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = length(init_params))
  theta_vals[1, ] <- init_params

  b_vals <- matrix(NA, nrow = iters+1, ncol = length(init_bs))
  b_vals[1, ] <- init_bs

  acceptance <- rep(NA, iters)

  log_post <- log_posterior(init_params, y, b_vals[1, ], stress, runout)

  # Progress bar in console
  pb <- txtProgressBar(min = 0, max = iters, style = 3)


  # MCMC loop
  for (k in seq_len(iters)) {


    prop <- MASS::mvrnorm(1, mu = c(theta_vals[k, ], b_vals[k, ]), Sigma = (tuner^2 * cov_matrix))
    theta_prop <- prop[1:length(init_params)]
    b_prop <- prop[(length(init_params)+1):length(prop)]


    log_post_prop <- log_posterior(theta_prop, y, b_prop, stress, runout)

    accept_prob <- exp(log_post_prop - log_post)

    if (accept_prob > runif(1)) {

      theta_vals[k+1, ] <- theta_prop
      b_vals[k+1, ] <- b_prop
      log_post <- log_post_prop
      acceptance[k] <- TRUE

    } else {

      theta_vals[k+1, ] <- theta_vals[k, ]
      b_vals[k+1, ] <- b_vals[k, ]
      acceptance[k] <- FALSE
    }


    setTxtProgressBar(pb, value = k)
  }


  if (show_plot) {
    par(mfrow = c(2, 3))
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

  list(theta = theta_vals, b = b_vals)
}



zz <- mcmc_mh_cov_f(5000, 1000, tail(pilot$theta, 1), tail(pilot$b, 1), cov_D, 0.35, fatigue$N, fatigue$s, fatigue$ro)
