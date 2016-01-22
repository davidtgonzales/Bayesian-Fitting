source('./process_data.R')
source('./plotting.R')
source('./gompertz.R')
setwd('../data/')

plotall<-function(startvalues,limits,triplicate,timerange,OD,times,names){
  layout(matrix(c(1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9), 4, 4, byrow = TRUE))
  par(mar=c(4.5,5,2,1))
  mlecol='red'
  bayescol='dodgerblue3'
  PlotCurve(log(OD[triplicate[1]:triplicate[3],timerange[1]:timerange[2]]),times[timerange[1]:timerange[2]],'grey',limits,'yes','time(minutes)','ln(OD540)','DH5A')
  growth<-c()
  value<-c(t(log(OD[triplicate[1],timerange[1]:timerange[2]])),t(log(OD[triplicate[2],timerange[1]:timerange[2]])),t(log(OD[triplicate[3],timerange[1]:timerange[2]])))
  level<-rep(c(times)[timerange[1]:timerange[2]],3)
  growth<-data.frame(value,level)
  colnames(growth)<-c('value','level')
  dev<-tapply(growth$value,growth$level,sd)
  growth$dev<-rep(dev,3)
  mlefit<-MaxLikelihood(growth,startvalues)
  chain<-run_metropolis_MCMC(mlefit,100000,growth$level,growth$value,growth$dev)
  burnIn=10000
  mcmcfit<-c(median(chain[-burnIn:-1,1]),median(chain[-burnIn:-1,2]),median(chain[-burnIn:-1,3]),median(chain[-burnIn:-1,4]))
  x<-seq(min(growth$level),max(growth$level),0.1)
  lines(x,Gompertz(mlefit,x),col=mlecol,lty=1,lwd=1)
  mq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.5))})
  uq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.975))})
  lq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.025))})
  lines(x,Gompertz(mq,x),col=bayescol,lty=1,lwd=1)
  lines(x,Gompertz(lq,x),col=bayescol,lty=2,lwd=1)
  lines(x,Gompertz(uq,x),col=bayescol,lty=2,lwd=1)
  legend('topleft',c('MLE','Bayesian'),col=c('red','dodgerblue3'),lty=c(1,1),box.lty=0,bty='n',lwd=2,cex=1.5)
  par(mar=c(2,5,2,1))
  hist(chain[-burnIn:-1,1],nclass=30,main="mu",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,2],nclass=30,main="K",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,3],nclass=30,main="lambda",col='dodgerblue3',border=NA,ylab='',xlab='')
  hist(chain[-burnIn:-1,4],nclass=30,main="lnN0",col='dodgerblue3',border=NA,ylab='',xlab='')
  par(mar=c(3,5,1,1))
  plot(chain[,1],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,2],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,3],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  plot(chain[,4],type="l",main="",col='dodgerblue3',ylab='',xlab='')
  print(mcmcfit)
}


OD<-read.csv('./dh5a_abs.csv',sep=',',header=TRUE)[2:101]
times<-t(read.csv('./dh5a_abs.csv',sep=',',header=FALSE)[1,2:101])
names<-read.csv('./dh5a_abs.csv',sep=',',header=TRUE)[1]

#DH5A uninduced
startvalues<-c(0.001,0.65,180,-0.65)
limits<-c(0,1000,-0.75,0.05)
triplicate<-c(1,2,3)
timerange<-c(1,100)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#DH5A induced
startvalues<-c(0.001,0.65,180,-0.65)
limits<-c(0,1000,-0.75,0.05)
triplicate<-c(4,5,6)
timerange<-c(1,100)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

#PatoG-OG241 in DH5A uninduced
startvalues<-c(0.001,0.5,60,-0.60)
limits<-c(0,600,-0.75,0.02)
triplicate<-c(7,8,9)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoG-OG241 in DH5A induced
startvalues<-c(0.001,0.5,60,-0.60) 
limits<-c(0,600,-0.75,0.02)
triplicate<-c(10,11,12)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

#PatoL-OG241 in DH5A uninduced
startvalues<-c(0.002,0.5,50,-0.15)
limits<-c(0,600,-0.2,0.5)
triplicate<-c(25,26,27)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoL-OG241 in DH5A induced
startvalues<-c(0.002,0.5,50,-0.15)
limits<-c(0,600,-0.2,0.5)
triplicate<-c(46,47,48)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

#PatoM-OG241 in DH5A uninduced
startvalues<-c(0.0016,0.52,75,-0.15)
limits<-c(0,600,-0.2,0.5)
triplicate<-c(49,50,51)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoM-OG241 in DH5A induced
startvalues<-c(0.0016,0.52,75,-0.15)
limits<-c(0,600,-0.2,0.5)
triplicate<-c(70,71,72)
timerange<-c(1,60)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

OD<-read.csv('./dh5a_nissle.csv',sep=',',header=TRUE)[102:201]
times<-t(read.csv('./dh5a_nissle.csv',sep=',',header=FALSE)[1,102:201])
names<-read.csv('./dh5a_nissle.csv',sep=',',header=TRUE)[1]

#PatoL-LOG241 in DH5A uninduced
startvalues<-c(0.001,0.65,1,-0.65)
limits<-c(0,1000,-0.7,0)
triplicate<-c(7,8,9)
timerange<-c(1,100)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoL-LOG241 in DH5A induced
startvalues<-c(0.001,0.65,120,-0.65)
limits<-c(0,1000,-0.7,0)
triplicate<-c(10,11,12)
timerange<-c(1,100)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

#PatoL-OG241 in Nissle uninduced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(13,14,15)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoL-OG241 in Nissle induced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(16,17,18)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

#PatoL-LOG241 in Nissle uninduced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(19,20,21)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#PatoL-LOG241 in Nissle induced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(22,23,24)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)

OD<-read.csv('./nissle_abs.csv',sep=',',header=TRUE)[2:101]
times<-t(read.csv('./nissle_abs.csv',sep=',',header=FALSE)[1,2:101])
names<-read.csv('./nissle_abs.csv',sep=',',header=TRUE)[1]

#Nissle uninduced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(31,32,33)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)
#Nissle induced
startvalues<-c(0.005,0.50,60,-0.60)
limits<-c(0,300,-0.75,0.1)
triplicate<-c(34,35,36)
timerange<-c(1,30)
plotall(startvalues,limits,triplicate,timerange,OD,times,names)