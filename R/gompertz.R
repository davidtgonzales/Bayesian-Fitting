Gompertz<-function(parameters,t){
  mu_m<-parameters[1]
  A<-parameters[2]
  lambda<-parameters[3]
  lnNo<-parameters[4]
  y<-lnNo+A*exp(-exp((((mu_m*exp(1))/A)*(lambda-t))+1))
  return(y)
}

LogLikelihood<-function(parameters,growth){
  #Gets the sum of the log likelihoods from a given set of parameters against the data. 
  logl<-0
  rows<-nrow(growth)
  pred<-Gompertz(parameters,growth$level)
  for(i in 1:rows){
    logl<-logl+dnorm(growth$value[i],mean=pred[i],sd=growth$dev[i],log=TRUE)
  }
  return(logl)
}

MaxLikelihood<-function(growth,parameters){
  #Gets the parameters from MLE.
  FitLikelihood<-function(growth,parameters){
    f<-function(w){-LogLikelihood(w,growth)}
    v<-optim(par=parameters,fn=f,control=list(maxit=1000))
    return(v)
  }
  mle<-FitLikelihood(growth,parameters)
  if(mle$convergence!= 0){
    cat("warning: no convergence in mle")
  }   
  return(mle$par)
}

#Likelihood function
likelihood<-function(param,x,y,dev){
  pred=Gompertz(param,x)
  singlelikelihoods=dnorm(y,mean=pred,sd=dev,log=T)
  sumll=sum(singlelikelihoods)
  return(sumll)
}

# Prior distribution
prior <- function(param){
  mu_m<-param[1]
  A<-param[2]
  lambda<-param[3]
  lnNo<-param[4]
  mu_mprior = dunif(mu_m, min=0, max=0.1, log = T)
  Aprior = dunif(A, min=0, max=1, log = T)
  lambdaprior = dunif(lambda, min=0, max=200, log = T)
  lnNoprior = dunif(lnNo, min=-1, max=0, log = T)
  return(lnNoprior+Aprior+mu_mprior+lambdaprior)
}

posterior <- function(param,x,y,dev){
  return (likelihood(param,x,y,dev) + prior(param))
}

proposalfunction <- function(param){
  return(rnorm(4,mean=param,sd=c(0.0001,0.01,1,0.01)))
}

run_metropolis_MCMC<-function(startvalue,iterations,x,y,dev){
  chain=array(dim=c(iterations+1,4))
  chain[1,]=startvalue
  for (i in 1:iterations){
    proposal=proposalfunction(chain[i,])
    probab=exp(posterior(proposal,x,y,dev)-posterior(chain[i,],x,y,dev))
    if (any(is.na(probab))==FALSE){
      if (runif(1) < probab){
        chain[i+1,]=proposal
      }else{
        chain[i+1,]=chain[i,]
      }
    }else{
      chain[i+1,]=chain[i,]
    }
  }
  return(chain)
}
