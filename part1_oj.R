
set.seed(101010)

# Import data
rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

# Try to see whether there's any immediately obvious effect from treatment
library(ggplot2)
ggplot(rats, aes(group = factor(rx), fill = factor(rx), x = time)) + geom_density(alpha = 0.3)
ggplot(rats, aes(x = factor(rx), y = time, fill = factor(rx))) + geom_violin()


# Negative log-likelihood using values from built-in Weibull functions
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








# Set up model matrices
rats$litter <- factor(rats$litter)

X <- model.matrix(~ 1 + rx, data = rats)
Z <- model.matrix(~ litter - 1, data = rats)
y <- rats[, c("time", "status")]



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






# Log-likelihood function l(y) = l(y, b)
lfyb <- function(theta, y, b, X, Z) {

  beta <- theta[1:2]
  log_sigma <- theta[3]
  log_sigma_b <- theta[4]

  eta <- X%*%beta + Z%*%b

  status <- y[, "status"]
  time <- y[, "time"]

  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_re_prob(time, eta, log_sigma),
    (1-status) * weib_re_surv(time, eta, log_sigma)
  )

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb


  # Gradient of log-likelihood with respect to b
  # 2 terms: d(lfy_b)/db + d(lfb)/db
  g <- t(Z) %*% (status*attr(weib_re_prob(time, eta, log_sigma), "gradient") + (1-status)*attr(weib_re_surv(time, eta, log_sigma), "gradient")) - b/(exp(log_sigma_b)^2)

  # Hessian of log-likelihood wrt b i.e. dl/(db db^T)) will be diagonal since
  # taking derivative wrt bi then wrt bj where j!=i gives us 0
  # i.e. we can use "hessian" from weib_me_*()
  # So H_diag is the diagonal entries from the Hessian matrix (and all other entries are 0)
  H_diag <- t(Z) %*% (status*attr(weib_re_prob(time, eta, log_sigma), "hessian") + (1-status)*attr(weib_re_surv(time, eta, log_sigma), "hessian")) - 1/(exp(log_sigma_b)^2)

  list(lf = lf, g = g, H_diag = H_diag)

}



# Check that output seems OK
lfyb(rep(1, 4), y, rep(0, ncol(Z)), X, Z)




# Laplace approximation calculation of marginal negative log-likelihood of theta
lal <- function(theta, y, X, Z) {

  # Starting guess for b
  b_init <- 0.01*rnorm(ncol(Z))


  opt <- optim(
    b_init,
    fn = function(b) {
      tmp <- lfyb(theta, y, b, X, Z)
      tmp$lf
    },
    gr = function(b) {
      tmp <- lfyb(theta, y, b, X, Z)
      tmp$g
    },
    method = "BFGS",

    # Set `fnscale = -1` so that we MINIMISE the NEGATIVE log-likelihood
    control = list(fnscale = -1)
  )


  b_opt <- opt$par

  soln_opt <- lfyb(theta, y, b_opt, X, Z)

  # Hessian was diagonal, so det(-H) is product of diagonal elements
  # i.e. log(det(-H)) is sum of logs of diagonal elements
  # We use abs() since sometimes very small values have the wrong sign due to floating point errors
  log_det_H <- sum(log(abs(soln_opt$H_diag)))


  # Return approximation of log-likelihood:
  # (Use Laplace approximation to find approx marginal likelihood (integral of
  # f(y, b) over b), then take the log; we do all that in one go)
  lap <- (length(b_opt)/2)*log(2*pi) + soln_opt$lf - log_det_H/2


  # Return negative log-likelihood, and also the optimal b that gives this
  attr(lap, "b") <- b_opt
  -lap
}


# Check that output is OK
lal(rep(1, 4), y, X, Z)



# Now we want to minimise negative log-likelihood (as returned by lal())
# Guess some starting values
# TODO justfy these
theta_init <- c("beta0" = 5,
                "beta1" = -1,
                "log_sigma" = -1,
                "log_sigma_b" = -2)


opt <- optim(
  theta_init,
  fn = lal,
  y = rats[, c("time", "status")],
  X = model.matrix(~ 1 + rx, data = rats),
  Z = model.matrix(~ litter - 1, data = rats),
  method = "BFGS",
  hessian = TRUE,
  control = list(trace = 1, maxit = 1000)
)



# Calculate standard errors from inverse Hessian
std_err <- sqrt(diag(solve(opt$hessian)))

# 95% confidence interval for each parameter
# Note asymptotic distribution is normal
data.frame(
  val = opt$par,
  se = std_err,
  lower = opt$par - qnorm(0.975)*std_err,
  upper = opt$par + qnorm(0.975)*std_err
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





mcmc_mh <- function(iters, burnin, init_params, init_bs, tuners, b_tuner, y, X, Z, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = length(init_params))
  colnames(theta_vals) <- names(init_params)
  theta_vals[1, ] <- init_params

  b_vals <- matrix(NA, nrow = iters+1, ncol = length(init_bs))
  b_vals[1, ] <- init_bs

  acceptance <- rep(NA, iters)
  acceptance_b <- rep(NA, iters)

  log_post <- log_posterior(init_params, y, init_bs, X, Z)

  # Progress bar in console
  pb <- txtProgressBar(min = 0, max = iters, style = 3)


  # MCMC loop
  for (k in seq_len(iters)) {


    # "Non-random" parameters
    theta_prop <- rnorm(length(init_params), mean = theta_vals[k, ], sd = tuners)


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
    b_prop <- rnorm(length(init_bs), mean = b_vals[k, ], sd = b_tuner)

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

    rows <- floor(sqrt(ncol(theta_vals)))
    par(mfrow = c(rows, rows+(ncol(theta_vals)%%2)))

    lapply(names(init_params), function(nm) {
      y <- theta_vals[, nm]
      plot(y, type = "l", main = nm, xlab = "Iteration", ylab = nm)
      rect(0, min(y), burnin, max(y), density = 10, col = "red")
    })

    par(mfrow = c(1, 1))
  }

  close(pb)
  cat(sprintf("Total acceptance:       %2.3f%%\n", mean(acceptance)*100))
  cat(sprintf("Burned-in acceptance:   %2.3f%%\n", mean(acceptance[-(1:burnin)])*100))
  cat(sprintf("Burned-in b acceptance: %2.3f%%\n", mean(acceptance_b[-(1:burnin)])*100))

  list(theta = theta_vals, b = b_vals)
}




mcmc_mh_cov <- function(iters, burnin, init_params, init_bs, cov_matrix, tuner, y, X, Z, show_plot = TRUE) {

  theta_vals <- matrix(NA, nrow = iters+1, ncol = length(init_params))
  colnames(theta_vals) <- names(init_params)
  theta_vals[1, ] <- init_params

  b_vals <- matrix(NA, nrow = iters+1, ncol = length(init_bs))
  b_vals[1, ] <- init_bs

  acceptance <- rep(NA, iters)

  log_post <- log_posterior(init_params, y, b_vals[1, ], X, Z)

  # Propositions for all iterations
  props <- MASS::mvrnorm(iters, mu = c(init_params, init_bs), Sigma = tuner^2 * cov_matrix)


  # Progress bar in console
  pb <- txtProgressBar(min = 0, max = iters, style = 3)


  # MCMC loop
  for (k in seq_len(iters)) {

    # Use proposed values for theta and b
    theta_prop <- props[k, 1:length(init_params)]
    b_prop <- props[k, (length(init_params)+1):ncol(props)]

    log_post_prop <- log_posterior(theta_prop, y, b_prop, X, Z)

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

    rows <- floor(sqrt(ncol(theta_vals)))
    par(mfrow = c(rows, rows+(ncol(theta_vals)%%2)))

    lapply(names(init_params), function(nm) {
      y <- theta_vals[, nm]
      plot(y, type = "l", main = nm, xlab = "Iteration", ylab = nm)
      rect(0, min(y), burnin, max(y), density = 10, col = "red")
    })

    par(mfrow = c(1, 1))
  }

  close(pb)
  cat(sprintf("Total acceptance:       %2.3f%%\n", mean(acceptance)*100))
  cat(sprintf("Burned-in acceptance:   %2.3f%%\n", mean(acceptance[-(1:burnin)])*100))

  list(theta = theta_vals, b = b_vals)
}



iters <- 100000
burnin <- 2000
pilot <- mcmc_mh(
  iters, burnin,
  c(beta0 = 4, beta1 = 0, log_sigma = 0, log_sigma_b = -1),
  rep(0, 50), c(0.1, 0.1, 0.1, 0.1), 0.03,
  rats[, c("time", "status")], X, Z
)


D <- cbind(pilot$theta, pilot$b)[-(1:burnin), ]
psych::pairs.panels(tail(D, 5000)[, 1:ncol(pilot$theta)], pch = ".")


cov_D <- cov(D)

iters <- 50000
burnin <- 1000

adjusted <- mcmc_mh_cov(
  iters, burnin,
  drop(tail(pilot$theta, 1)),
  drop(tail(pilot$b, 1)),
  cov_D, 0.1,
  rats[, c("time", "status")], X, Z
)


adjD <- cbind(adjusted$theta, adjusted$b)[-(1:burnin), ]
psych::pairs.panels(tail(adjD, 5000)[, 1:ncol(adjusted$theta)], pch = ".")


par(mfrow = c(2, 2))

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
