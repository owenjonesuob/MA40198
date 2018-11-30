rats <- read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

nll <- function(theta){
  
  l1 = dweibull(rats$time[rats$status==1&rats$rx==1], shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]), log=TRUE) 
  l2 = dweibull(rats$time[rats$status==1&rats$rx==0], shape=1/exp(theta[3]), scale=exp(theta[1]), log=TRUE) 
  
  l3 = pweibull(rats$time[rats$status==0&rats$rx==1], shape=1/exp(theta[3]), scale=exp(theta[1] + theta[2]), log.p=TRUE, lower.tail=FALSE) 
  l4 = pweibull(rats$time[rats$status==0&rats$rx==0], shape=1/exp(theta[3]), scale=exp(theta[1]), log.p=TRUE, lower.tail=FALSE)   
  
  -(sum(l1)+sum(l2)+sum(l3)+sum(l4))

}


theta0=c(1,1,1)

q=optim(par=theta0,fn=nll, method="BFGS", hessian=TRUE)
inv_hess = solve(q$hessian)

#2. 
q$par[2] + 1.96*sqrt(inv_hess[2,2])
q$par[2] - 1.96*sqrt(inv_hess[2,2])

#Interpretation: This shows that at the 5% significance level, beta1 is significantly different from 0 and therefore 
#we can say that there is sufficient evidence to suggest that recieving treatment has an effect on the time it takes
#for a tumour to appear. We have a negative estimate, which decreases the scale parameter in the Weibull distribution. 
#This in turn DECREASES/INCREASES the esimated times for tumour appearance. 

#3. 

install.packages("actuar")
library(actuar)

nll_ll <- function(theta){
  
  shape = 1/exp(theta[3])
  scale_0 = exp(theta[1])
  scale_1 = exp(theta[1]+theta[2])
  
  l1 = log(shape)-log(scale_1)+(shape-1)*log(rats$time[rats$status==1&rats$rx==1]/scale_1) - 2*log(1+(rats$time[rats$status==1&rats$rx==1]/scale_1)^shape)
  l2 = log(shape)-log(scale_0)+(shape-1)*log(rats$time[rats$status==1&rats$rx==0]/scale_0) - 2*log(1+(rats$time[rats$status==1&rats$rx==0]/scale_0)^shape)
  
  l3 = -log(1+(rats$time[rats$status==0&rats$rx==1]/scale_1)^shape)
  l4 = -log(1+(rats$time[rats$status==0&rats$rx==0]/scale_0)^shape)
  
  -(sum(l1)+sum(l2)+sum(l3)+sum(l4))
  
}


theta0_ll=c(1,1,1)
q_ll=optim(par=theta0_ll,fn=nll_ll, method="BFGS", hessian=TRUE)
inv_hess_ll = solve(q_ll$hessian)
q_ll$par[2] + 1.96*sqrt(inv_hess_ll[2,2])
q_ll$par[2] - 1.96*sqrt(inv_hess_ll[2,2])


#4.
rats$litter <- factor(rats$litter)
Z <- model.matrix(~ litter - 1,rats)

X <- model.matrix(~ rx, rats)


lfyb <- function(theta, y, b, X, Z){
  eta_0 = as.numeric(theta[1] + Z%*%b)
  eta_1 = as.numeric(X%*%theta[1:2] + Z%*%b)
  scale_0 = exp(eta_0)
  scale_1 = exp(eta_1)
  shape = 1/exp(theta[3])
  
  shape = 1/exp(theta[3])

  l1 = dweibull(y[rats$status==1&rats$rx==1], shape=shape, scale=scale_1, log=TRUE) 
  l2 = dweibull(y[rats$status==1&rats$rx==0], shape=shape, scale=scale_0, log=TRUE) 
  
  l3 = pweibull(y[rats$status==0&rats$rx==1], shape=shape, scale=scale_1, log.p=TRUE, lower.tail=FALSE) 
  l4 = pweibull(y[rats$status==0&rats$rx==0], shape=shape, scale=scale_0, log.p=TRUE, lower.tail=FALSE)   
  
  lfy_b = sum(l1)+sum(l2)+sum(l3)+sum(l4)
  lfb = sum(dnorm(b, 0, exp(theta[4]), log=TRUE))
  -lfy_b - lfb
  
  
}

lal <- function(theta, y, X, Z){
  
  b <- rep(0,ncol(Z))
  opt = optim(par=b, lfyb, theta=theta, y=y, X=X, Z=Z, method="BFGS", hessian=TRUE)
  
  
  la <- -lfyb(theta, y, opt$par, X, Z) + length(b)*log(2*pi)/2 - sum(log(abs(diag(solve(opt$hessian)))))/2
  
  attr(la,"b") <- as.numeric(b)
  
  -la
  
}

theta0=c(1.1749184,  0.6614587, -1.4266831,  1.7078852)

optim(par=theta0, lal, X=X, y=rats$time, Z=Z, method="BFGS", hessian=TRUE)

#works w/ nelder mead
