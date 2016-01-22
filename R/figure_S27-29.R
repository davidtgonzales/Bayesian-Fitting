source('./process_data.R')
source('./plotting.R')
source('./hill.R')
setwd('../data/induction')
address=getwd()

groups<-c('patog','patol','patom')
factor='acetoacetate'
data<-ProcessData(address)
medians<-GetMedians(data)

for (group in groups){
  mediansubset<-GetMedianSubset(medians,group,factor) 
  #MLE & MCMC fitting
  if (group=='patog'){
    startvalues<-c(700,12000,0.15,2.5)
  }
  if (group=='patol'){
    startvalues<-c(600,4500,0.2,1.5)
  }
  if (group=='patom'){
    startvalues<-c(600,3000,0.2,1.5)
  }
  mlefit<-MaxLikelihood(mediansubset,startvalues)
  chain<-run_metropolis_MCMC(mlefit,100000,mediansubset$level,mediansubset$value,mediansubset$dev)
  #Plot results
  burnIn=10000
  mcmcfit<-c(median(chain[-burnIn:-1,1]),median(chain[-burnIn:-1,2]),median(chain[-burnIn:-1,3]),median(chain[-burnIn:-1,4]))
  print(group)
  print(mcmcfit)
  layout(matrix(c(1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9), 4, 4, byrow = TRUE))
  par(mar=c(4.5,5,2,1))
  plot(mediansubset$level,mediansubset$value,log='x',ylab='RFU',xlab='mM aceoacetate',type='p',pch=20,col='grey',cex.axis=1.5,cex.lab=1.7,cex=3)
  grid(NULL, NULL, lwd = 1)
  x<-seq(min(mediansubset$level),max(mediansubset$level),0.0001)
  lines(x,Hill(mlefit,x),col='red',lty=1,lwd=1)
  mq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.5))})
  uq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.975))})
  lq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.025))})
  lines(x,Hill(mq,x),col='dodgerblue3',lty=1,lwd=1)
  lines(x,Hill(lq,x),col='dodgerblue3',lty=2,lwd=1)
  lines(x,Hill(uq,x),col='dodgerblue3',lty=2,lwd=1)
  legend('topleft',c('MLE','Bayesian'),col=c('red','dodgerblue3'),lty=c(1,1),box.lty=0,bty='n',lwd=2,cex=1.5)
  par(mar=c(2,4,2,1))
  hist(chain[-burnIn:-1,1],nclass=30,main="Pmin",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,2],nclass=30,main="Pmax",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,3],nclass=30,main="Kd",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,4],nclass=30,main="n",col='dodgerblue3',border=NA,ylab='',xlab='')
  par(mar=c(3,4,2,1))
  plot(chain[,1],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,2],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,3],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,4],type="l",main="",col='dodgerblue3',ylab='',xlab='')
}
