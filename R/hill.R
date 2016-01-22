Hill<-function(parameters,x){
  #Calculates y of the hill function given parameters and x values.
  pmin<-parameters[1]
  pmax<-parameters[2]
  Kd<-parameters[3]
  n<-parameters[4]
  y<-pmin+(pmax-pmin)*((x/Kd)^n)/(1+(x/Kd)^n)
  return(y)
}

LogLikelihood<-function(parameters,mediansubset){
  #Gets the sum of the log likelihoods from a given set of parameters against the data. 
  logl<-0
  rows<-nrow(mediansubset)
  pred<-Hill(parameters,mediansubset$level)
  for(i in 1:rows){
    logl<-logl+dnorm(mediansubset$value[i],mean=pred[i],sd=mediansubset$dev[i],log=TRUE)
  }
  return(logl)
}

MaxLikelihood<-function(mediansubset,parameters){
  #Gets the parameters from MLE.
  FitLikelihood<-function(mediansubset,parameters){
    f<-function(w){-LogLikelihood(w,mediansubset)}
    v<-optim(par=parameters,fn=f,control=list(maxit=1000))
    return(v)
  }
  mle<-FitLikelihood(mediansubset,parameters)
  if(mle$convergence!= 0){
    cat("warning: no convergence in mle")
  }
  return(mle$par)
}

#Log likelihood function
likelihood<-function(param,x,y,dev){
  pred=Hill(param,x)
  singlelikelihoods=dnorm(y,mean=pred,sd=dev,log=T)
  sumll=sum(singlelikelihoods)
  return(sumll)
}

# Prior distribution
prior <- function(param){
  pmin<-param[1]
  pmax<-param[2]
  Kd<-param[3]
  n<-param[4]
  pminprior = dunif(pmin, min=0, max=1000, log = T)
  pmaxprior = dunif(pmax, min=0, max=20000, log = T)
  Kdprior = dunif(Kd, min=0, max=1, log = T)
  nprior = dunif(n, min=0, max=4, log = T)
  return(pminprior+pmaxprior+Kdprior+nprior)
}

posterior <- function(param,x,y,dev){
  return (likelihood(param,x,y,dev) + prior(param))
}

proposalfunction <- function(param){
  return(rnorm(4,mean=param,sd=c(1,10,0.01,0.1)))
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