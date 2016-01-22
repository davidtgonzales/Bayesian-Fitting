source('./process_data.R')
source('./plotting.R')
source('./hill.R')
setwd('../data/induction')
address=getwd()
data<-ProcessData(address)
#PatoG FACS plot
group='patog'
factor='acetoacetate'
medians<-GetMedians(data)
mediansubset<-GetMedianSubset(medians,group,factor) 
#colors<-rev(heat.colors(8))
#library(RColorBrewer)
#colors<-brewer.pal(8,'YlGnBu')
colors<-c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
par(mfrow=c(2,2))
par(mar=c(4.5,5,4,1))
PlotHist(data,group,'acetoacetate','mM',c(2,5,0,2.5),'none',colors)
par(adj=0)
title('A)',cex.main=1.7)
par(adj=0.5)
i=1
y<-c(2.6,2.5,2.4,2.3,2.2,2.1,2.0,1.9)
meanpoints<-c()
for (level in unique(medians$level)){
  mid<-mean(log10(medians[medians$group==group & medians$level==level,]$median))
  stdev<-sd(log10(medians[medians$group==group & medians$level==level,]$median))
  meanpoints<-c(meanpoints,mid)
  segments(mid-stdev,y[i],mid+stdev,y[i],col=colors[i])
  epsilon<-0.1
  segments(mid-stdev,y[i]-epsilon,mid-stdev,y[i]+epsilon,col=colors[i])
  segments(mid+stdev,y[i]-epsilon,mid+stdev,y[i]+epsilon,col=colors[i])
  i=i+1
}
points(meanpoints,c(2.6,2.5,2.4,2.3,2.2,2.1,2.0,1.9),col=colors,pch=20,cex=1)
#PatoG characterisation curve
x<-mediansubset$level
y<-mediansubset$value
plot(x,y, main='',xlab='mM Acetoacetate',ylab='RFU',type='p',pch=20,col=rep(colors,each=3),log='x',cex=2,ylim=c(0,20000),cex.lab=1.7,cex.axis=1.5)
grid(NULL, NULL, lwd = 1)
startvalues<-c(700,12000,0.15,2.5) #patog
mlefit<-MaxLikelihood(mediansubset,startvalues)
chain<-run_metropolis_MCMC(mlefit,100000,mediansubset$level,mediansubset$value,mediansubset$dev)
burnIn=10000
mcmcfit<-c(median(chain[-burnIn:-1,1]),median(chain[-burnIn:-1,2]),median(chain[-burnIn:-1,3]),median(chain[-burnIn:-1,4]))
x<-seq(min(mediansubset$level),max(mediansubset$level),0.0001)
mq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.5))})
uq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.975))})
lq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.025))})
lines(x,Hill(mq,x),col='black',lty=1,lwd=1)
lines(x,Hill(lq,x),col='black',lty=2,lwd=1)
lines(x,Hill(uq,x),col='black',lty=2,lwd=1)
legend('topleft',c('0mM acetoacetate','1e-3mM','1e-2mM','1e-1.5mM','1e-1mM','1e-0.5mM','1mM','10mM'),col=colors,box.lty=0,bty='n',pch=20,cex=0.8)
par(adj=0)
title('B)',cex.main=1.7)
par(adj=0.5)
#PatoL FACS plot
group='patol'
factor='acetoacetate'
medians<-GetMedians(data)
mediansubset<-GetMedianSubset(medians,group,factor) 
PlotHist(data,group,'acetoacetate','mM',c(2,5,0,2.5),'none',colors)
par(adj=0)
title('C)',cex.main=1.7)
par(adj=0.5)
i=1
y<-c(2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8)
meanpoints<-c()
for (level in unique(medians$level)){
  mid<-mean(log10(medians[medians$group==group & medians$level==level,]$median))
  stdev<-sd(log10(medians[medians$group==group & medians$level==level,]$median))
  meanpoints<-c(meanpoints,mid)
  segments(mid-stdev,y[i],mid+stdev,y[i],col=colors[i])
  epsilon<-0.1
  segments(mid-stdev,y[i]-epsilon,mid-stdev,y[i]+epsilon,col=colors[i])
  segments(mid+stdev,y[i]-epsilon,mid+stdev,y[i]+epsilon,col=colors[i])
  i=i+1
}
points(meanpoints,c(2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8),col=colors,pch=20,cex=1)
#PatoL characterisation curve
x<-mediansubset$level
y<-mediansubset$value
plot(x,y, main='',xlab='mM Acetoacetate',ylab='RFU',type='p',pch=20,col=rep(colors,each=3),log='x',cex=2,ylim=c(0,5000),cex.lab=1.7,cex.axis=1.5)
grid(NULL, NULL, lwd = 1)
startvalues<-c(600,4500,0.2,1.5) #patol and patom
mlefit<-MaxLikelihood(mediansubset,startvalues)
chain<-run_metropolis_MCMC(mlefit,100000,mediansubset$level,mediansubset$value,mediansubset$dev)
burnIn=10000
mcmcfit<-c(median(chain[-burnIn:-1,1]),median(chain[-burnIn:-1,2]),median(chain[-burnIn:-1,3]),median(chain[-burnIn:-1,4]))
x<-seq(min(mediansubset$level),max(mediansubset$level),0.0001)
mq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.5))})
uq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.975))})
lq<-apply(chain[-burnIn:-1,],MARGIN=2,function(x){quantile(x,p=c(0.025))})
lines(x,Hill(mq,x),col='black',lty=1,lwd=1)
lines(x,Hill(lq,x),col='black',lty=2,lwd=1)
lines(x,Hill(uq,x),col='black',lty=2,lwd=1)
#legend('topleft',c('0mM acetoacetate','1e-3mM','1e-2mM','1e-1.5mM','1e-1mM','1e-0.5mM','1mM','10mM'),col=colors,box.lty=0,bty='n',pch=20,cex=0.8)
par(adj=0)
title('D)',cex.main=1.7)
par(adj=0.5)

