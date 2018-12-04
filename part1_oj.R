
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

  status <- X[, 2]

  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_re_prob(y, eta, log_sigma),
    (1-status) * weib_re_surv(y, eta, log_sigma)
  )

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb


  # Gradient of log-likelihood with respect to b
  # 2 terms: d(lfy_b)/db + d(lfb)/db
  g <- t(Z) %*% (status*attr(weib_re_prob(y, eta, log_sigma), "gradient") + (1-status)*attr(weib_re_surv(y, eta, log_sigma), "gradient")) - b/(exp(log_sigma_b)^2)

  # Hessian of log-likelihood wrt b i.e. dl/(db db^T)) will be diagonal since
  # taking derivative wrt bi then wrt bj where j!=i gives us 0
  # i.e. we can use "hessian" from weib_me_*()
  # So H_diag is the diagonal entries from the Hessian matrix (and all other entries are 0)
  H_diag <- t(Z) %*% (status*attr(weib_re_prob(y, eta, log_sigma), "hessian") + (1-status)*attr(weib_re_surv(y, eta, log_sigma), "hessian")) - 1/(exp(log_sigma_b)^2)

  list(lf = lf, g = g, H_diag = H_diag)

}



# Check that output seems OK
lfyb(rep(1, 4), rats$time, rep(0, ncol(Z)), X, Z)




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
lal(rep(1, 4), rats$time, X, Z)





# Now we want to minimise negative log-likelihood (as returned by lal())
theta_init <- rep(1, 4)


# TODO attempt to find convergent starting value
max_guesses <- 1000

for (k in 1:max_guesses) {

  theta_init <- runif(4, min = -10, max = 10)

  tryCatch({

    opt <- optim(
      theta_init,
      fn = lal,
      gr = function(theta, y, X, Z) {
        grad(lal, theta, y = y, X = X, Z = Z)
      },
      y = rats$time,
      X = model.matrix(~ 1 + rx, data = rats),
      Z = model.matrix(~ litter - 1, data = rats),
      method = "BFGS",
      control = list(maxit = 100)
    )

    beepr::beep("fanfare")
    break

  },
  error = function(e) cat(k, "| Failed to converge: ", theta_init, "\n"))
}; beepr::beep("mario")











## Attempt to derive log likelihood gradient wrt theta...

# Expressions for model including random effects
weib_re_prob <- deriv(
  expression(
    log((1/exp(log_sigma))/exp(beta0 + beta1*status + b) *
          (t/exp(beta0 + beta1*status + b))^((1/exp(log_sigma))-1) *
          exp(-(t/exp(beta0 + beta1*status + b))^(1/exp(log_sigma))))
  ),
  namevec = c("beta0", "beta1", "log_sigma", "b"),
  function.arg = c("t", "beta0", "beta1", "status", "b", "log_sigma"),
  hessian = TRUE
)

weib_re_surv <- deriv(
  expression(
    -(t/exp(beta0 + beta1*status + b))^(1/exp(log_sigma))
  ),
  namevec = c("beta0", "beta1", "log_sigma", "b"),
  function.arg = c("t", "beta0", "beta1", "status", "b", "log_sigma"),
  hessian = TRUE
)



# Log-likelihood function l(y) = l(y, b)
lfyb <- function(theta, y, b, X, Z) {

  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]
  log_sigma_b <- theta[4]

  b_ <- Z%*%b

  status <- X[, "rx"]

  # Log conditional density of y given b
  lfy_b <- sum(
    status * weib_re_prob(y, beta0, beta1, status, b_, log_sigma),
    (1-status) * weib_re_surv(y, beta0, beta1, status, b_, log_sigma)
  )

  # Log marginal density of b
  lfb <- sum(dnorm(x = b, mean = 0, sd = exp(log_sigma_b), log = TRUE))

  # Log joint density of y and b is the sum (joint density is product - y, b are independent)
  lf <- lfy_b + lfb


  # Gradient of log-likelihood with respect to b
  # 2 terms: d(lfy_b)/db + d(lfb)/db
  g <- t(Z) %*% (status*attr(weib_re_prob(y, beta0, beta1, status, b_, log_sigma), "gradient")[, "b"] +
                   (1-status)*attr(weib_re_surv(y, beta0, beta1, status, b_, log_sigma), "gradient")[, "b"]) - b/(exp(log_sigma_b)^2)

  # Hessian of log-likelihood wrt b i.e. dl/(db db^T)) will be diagonal since
  # taking derivative wrt b_i then wrt b_j where j!=i gives us 0
  # i.e. we can use "hessian" from weib_me_*()
  # So H_diag is the diagonal entries from the Hessian matrix (and all other entries are 0)
  H_diag <- t(Z) %*% (status*attr(weib_re_prob(y, beta0, beta1, status, b_, log_sigma), "hessian")[, "b", "b"] + (1-status)*attr(weib_re_surv(y, beta0, beta1, status, b_, log_sigma), "hessian")[, "b", "b"]) - 1/(exp(log_sigma_b)^2)

  list(lf = lf, g = g, H_diag = H_diag)

}



# Check that output seems OK
lfyb(rep(1, 4), rats$time, rep(0, ncol(Z)), X, Z)




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
  # We use abs() to ENSURE elements are positive, because occasionally due to
  # floating point errors, we do get very small elements with the wrong sign
  log_det_H <- sum(log(abs(soln_opt$H_diag)))


  # Return approximation of log-likelihood:
  # (Use Laplace approximation to find approx marginal likelihood (integral of
  # f(y, b) over b), then take the log; we do all that in one go)
  lap <- (length(b_opt)/2)*log(2*pi) + soln_opt$lf - log_det_H/2





  beta0 <- theta[1]
  beta1 <- theta[2]
  log_sigma <- theta[3]
  log_sigma_b <- theta[4]
  status <- X[, "rx"]
  g <- colSums(
    t(Z) %*% (status*attr(weib_re_prob(y, beta0, beta1, status, b_opt, log_sigma), "gradient")[, c("beta0", "beta1", "log_sigma")] +
                (1-status)*attr(weib_re_surv(y, beta0, beta1, status, b_opt, log_sigma), "gradient")[, c("beta0", "beta1", "log_sigma")])
  )
  g <- c(g, "log_sigma_b" = sum(-b_opt/(log_sigma_b^2)))


  # Return negative log-likelihood, and also the optimal b that gives this
  attr(lap, "b") <- b_opt
  attr(lap, "gradient") <- g

  -lap
}


# Check that output is OK
lal(rep(1, 4), rats$time, X, Z)





# Now we want to minimise negative log-likelihood (as returned by lal())
theta_init <- rep(0, 4)

opt <- optim(
  theta_init,
  fn = lal,
  gr = function(x, y, X, Z) attr(lal(x, y, X, Z), "gradient"),
  y = rats$time,
  X = model.matrix(~ 1 + rx, data = rats),
  Z = model.matrix(~ litter - 1, data = rats),
  method = "BFGS",
  control = list(REPORT = 1, trace = 1, maxit = 100)
)





