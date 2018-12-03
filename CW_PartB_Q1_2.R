#2

fatigue<- read.table("http://people.bath.ac.uk/kai21/ASI/fatigue.txt")


#nll
nll <-function(theta){
  alpha<-exp(theta[1])
  delta<-theta[2]
  sigma<-exp(theta[3])
  gam<-80*(1/(1+exp(-theta[4])))
  
  
  #seperate datasets into 
  fatigue0=subset(fatigue,ro==0)
  fatigue1=subset(fatigue,ro==1)
  
  
  l1<-(dweibull(fatigue0$N,shape=(1/sigma),scale=(alpha*(fatigue0$s-gam)^(delta)),log =TRUE))
  
  l2<-(pweibull(fatigue1$N,shape=(1/sigma),lower.tail = FALSE,scale=(alpha*(fatigue1$s-gam)^(delta)),log.p =TRUE))
  
  
  nll<--sum(l1)-sum(l2)
  nll
}


theta=c(1,1,1,1)

nll(theta)

theta0=c(10,10,10,10)

p=optim(par=theta0,fn=nll,control=list(trace=1,REPORT=1),method="BFGS",hessian = TRUE)


alpha<-exp(p$par[1])
delta<-p$par[2]
sigma<-exp(p$par[3])
gam<-80*(1/(1+exp(-theta[4])))

#check
N<-(alpha*(fatigue$s-gam)^(delta))*rweibull(26,shape=1/sigma,scale=1)

plot(N)

plot(fatigue$N)


#confidence intervals

std_err=sqrt(diag(solve(p$hessian)))

names(std_err) <- c("log_alpha", "delta", "log_sigma")

p$par[1] + c(-1, 1)*qnorm(0.975)*std_err[2]
p$par[2] + c(-1, 1)*qnorm(0.975)*std_err[2]
p$par[3] + c(-1, 1)*qnorm(0.975)*std_err[2]
