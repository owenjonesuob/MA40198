#Script for question 1 and 2 on MA40198 Coursework

#Start by finding negative log likelihodd and grad for use in optim

rats<-read.table("http://people.bath.ac.uk/kai21/ASI/rats_data.txt")

#nll
nll <-function(theta){
  beta0<-theta[1]
  beta1<-theta[2]
  logsigma<-theta[3]
  
  logshape<--logsigma
  ratshape<-(1/exp(logsigma))
  
  #seperate datasets into 
  rats00<-subset(rats,rx<1 & status<1 )
  rats10<-subset(rats,rx>0 & status<1)
  rats01<-subset(rats,rx<1 & status>0)
  rats11<-subset(rats, rx >0 & status>0)
  
  l1<-(dweibull(rats01$time,shape=(1/exp(logsigma)),scale=(exp(beta0)),log =TRUE))
  l2<-(dweibull(rats11$time,shape=(1/exp(logsigma)),scale=(exp(beta0+beta1)),log =TRUE))
  
  l3<-(pweibull(rats00$time,shape=(1/exp(logsigma)),scale=(exp(beta0)),lower.tail = FALSE,log.p =TRUE))
  l4<-(pweibull(rats10$time,shape=(1/exp(logsigma)),scale=(exp(beta0+beta1)),lower.tail = FALSE,log.p =TRUE))

  nll<--sum(l1)-sum(l2)-sum(l3)-sum(l4)
  nll
}

theta0<-c(1,1,1)

q=optim(par=theta0,fn=nll,control=list(trace=1,REPORT=1),method="Nelder-Mead")

p=optim(par=q$par,fn=nll,control=list(trace=1,REPORT=1),method="BFGS",hessian=TRUE)


# Calculate standard errors from inverse Hessian

llog_std_err <- sqrt(diag(solve(p$hessian)))

